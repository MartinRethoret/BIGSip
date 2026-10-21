#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Command-line tool and interactive wizard to download data from the BIGSdb
platforms hosted by BIGSdb-Pasteur and PubMLST without writing REST API calls:

- isolate genomes (contigs FASTA), from the ``<organism>_isolates`` database;
- allele FASTA files of loci (a whole scheme or chosen loci), from the
  ``<organism>_seqdef`` database;
- scheme profiles (ST tables), for schemes that define a primary key.

The user only gives an organism name. BIGSip builds the database name and the
API URLs, checks them against the site catalogue, then runs each request
through BIGSdb_downloader, which handles the OAuth authentication required for
data produced after 1 January 2025:

    https://github.com/kjolley/BIGSdb_downloader

BIGSdb REST API reference:

    https://bigsdb.readthedocs.io/en/latest/rest.html

Requirements
------------
Python >= 3.9 (BIGSdb_downloader requirement). BIGSip only uses the standard
library; ``requests`` is optional (faster metadata reads and ``--public``
mode) and is already installed with ``bigsdb-downloader`` (through ``rauth``).

License: see the LICENSE file of the repository.
"""
from __future__ import annotations

import argparse
import csv
import datetime
import difflib
import getpass
import json
import logging
import os
import re
import shlex
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple, Union
from urllib.parse import quote

try:  # Optional: plain HTTP access to public metadata and --public mode.
    import requests
except ImportError:  # pragma: no cover
    requests = None

__version__ = "1.0.0"

log = logging.getLogger("bigsip")

# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------
SITES: Dict[str, str] = {
    "Pasteur": "https://bigsdb.pasteur.fr/api",
    "PubMLST": "https://rest.pubmlst.org",
}
SITE_WEB: Dict[str, str] = {
    "Pasteur": "https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl",
    "PubMLST": "https://pubmlst.org/bigsdb",
}
DEFAULT_SITE = "Pasteur"
DEFAULT_TIMEOUT = 600          # seconds; whole genomes can be slow to stream
METADATA_TIMEOUT = 60
DEFAULT_OUTDIR = "bigsdb_downloads"
DOWNLOAD_LOG = "download_log.tsv"

CONFIG_ENV = "BIGSIP_CONFIG"
USER_DIR = Path.home() / ".bigsip"
DEFAULT_CONFIG = USER_DIR / "config.json"
SHARED_ENV = "BIGSIP_SHARED"
SHARED_CONFIG_NAME = "bigsip_shared.json"
# Tokens are personal: {user} and {home} are expanded for the current user.
DEFAULT_TOKEN_DIR = "~/.bigsip/tokens"
SHARED_KEYS = ("downloader", "downloader_python", "default_site", "token_dir",
               "key_names", "shared_tokens")

POLICY_TEXT = (
    "Data remain subject to the usage policy of each site:\n"
    "  - Pasteur: https://bigsdb.pasteur.fr/policy/\n"
    "  - PubMLST: https://pubmlst.org/terms-conditions"
)
PROG = "python bigsip.py"


# --------------------------------------------------------------------------
# Domain model
# --------------------------------------------------------------------------
class DbKind(str, Enum):
    """The two BIGSdb database types used by BIGSip."""

    SEQDEF = "seqdef"        # sequence definitions: loci, alleles, schemes, profiles
    ISOLATES = "isolates"    # isolate records and their genome assemblies

    @property
    def label(self) -> str:
        return "genomes (isolates)" if self is DbKind.ISOLATES else "loci / profiles (seqdef)"

    @property
    def data(self) -> str:
        return "genomes" if self is DbKind.ISOLATES else "loci / profiles"


@dataclass
class Organism:
    """An organism of the site catalogue and its databases."""

    name: str
    description: str = ""
    seqdef: Optional[str] = None
    isolates: Optional[str] = None

    def db(self, kind: DbKind) -> Optional[str]:
        return self.seqdef if kind is DbKind.SEQDEF else self.isolates


@dataclass
class Scheme:
    """A BIGSdb scheme. ``has_profiles`` is None when it could not be determined."""

    id: str
    description: str = ""
    has_profiles: Optional[bool] = None


@dataclass
class SchemeInfo:
    id: str
    description: str
    loci: List[str]


@dataclass
class Job:
    """One file to download: ``key`` identifies it in reports (isolate ID, locus...)."""

    key: str
    url: str
    dest: Path
    fmt: str = "fasta"       # expected format, checked before the file is kept


@dataclass
class BatchResult:
    ok: List[str] = field(default_factory=list)
    skipped: List[str] = field(default_factory=list)
    failed: List[Tuple[str, str]] = field(default_factory=list)


class BigsdbError(RuntimeError):
    """Error while preparing or running a BIGSdb request."""


class FatalError(BigsdbError):
    """Error that would make every remaining request fail: stop the batch."""


class AuthError(FatalError):
    """Authentication problem (missing, revoked or rejected credentials)."""


# --------------------------------------------------------------------------
# bigsdb_downloader error messages -> actionable messages
# --------------------------------------------------------------------------
# Raised by the API when profiles are requested for a scheme without primary key.
NO_PK_HINT = ("this scheme has no primary key (no ST) so it has no profiles; "
              "its sequences are available with 'loci --scheme'")

# Patterns taken from bigsdb_downloader's main.py stderr messages.
AUTH_PATTERNS: Tuple[Tuple[str, str], ...] = (
    (r"Run interactively|No access token saved|No client credentials saved",
     "Your BIGSdb keys/tokens are missing or must be set up again for this site."),
    (r"New access token required|Invalid access token|verification",
     "Your access token was revoked or has expired: you need to authorise again."),
    (r"client is unauthorized",
     "Your client key (client id/secret) is rejected by the site (revoked or invalid)."),
    (r"token directory .* not writable|Failed to create token directory",
     "The token directory is not writable."),
)


def short_message(text: str, max_len: int = 200) -> str:
    """Return the ``message`` of a BIGSdb JSON error if present, else a truncated text."""
    text = (text or "").strip().replace("\n", " ")
    if re.search(r"does not have a primary key", text):
        return NO_PK_HINT
    match = re.search(r'"message"\s*:\s*"([^"]+)"', text)
    return match.group(1) if match else text[:max_len]


def explain_downloader_error(msg: str, site: str) -> Optional[FatalError]:
    """
    Map a raw bigsdb_downloader error to an exception the user can act on.

    Returns None for ordinary per-request errors (missing isolate, 404...), which
    are reported in the failure list instead of stopping the batch.
    """
    msg = msg or ""
    missing = re.search(r"No module named '([^']+)'", msg)
    if missing or msg.strip().startswith("Traceback"):
        detail = ("the Python module '{}' is missing".format(missing.group(1)) if missing
                  else "bigsdb_downloader crashed")
        return FatalError(
            "bigsdb_downloader does not work: {}.\n"
            "  -> Fix 1: install its dependencies in the Python being used "
            "(pip install bigsdb-downloader, or pip install rauth)\n"
            "  -> Fix 2: point to the 'bigsdb-downloader' executable of a venv "
            "(setup --downloader /path/venv/bin/bigsdb-downloader)\n"
            "  -> Fix 3: give the Python that has rauth (setup --python /path/venv/bin/python)\n"
            "  (end of message: {})".format(
                detail, short_message(msg.strip().splitlines()[-1], 150)))
    for pattern, text in AUTH_PATTERNS:
        if re.search(pattern, msg, flags=re.IGNORECASE):
            return AuthError(
                "{}\n  -> Run: {} auth --site {} --organism <organism>\n"
                "  (bigsdb_downloader message: {})".format(text, PROG, site, short_message(msg)))
    return None


# --------------------------------------------------------------------------
# Configuration: defaults < shared (admin) config < personal config
# --------------------------------------------------------------------------
Config = Dict[str, object]


def config_path(cli_value: Optional[str] = None) -> Path:
    """Path of the PERSONAL config of the current user."""
    if cli_value:
        return Path(cli_value).expanduser()
    if os.environ.get(CONFIG_ENV):
        return Path(os.environ[CONFIG_ENV]).expanduser()
    return DEFAULT_CONFIG


def shared_config_path() -> Path:
    """Path of the SHARED config, written once by the admin next to the script."""
    if os.environ.get(SHARED_ENV):
        return Path(os.environ[SHARED_ENV]).expanduser()
    return Path(__file__).resolve().parent / SHARED_CONFIG_NAME


def default_config() -> Config:
    return {
        "downloader": None,          # bigsdb_downloader .py script or executable
        "downloader_python": None,   # interpreter for a .py downloader (needs rauth)
        "token_dir": None,           # personal directory for keys and OAuth tokens
        "default_site": DEFAULT_SITE,
        "key_names": {site: site for site in SITES},
        "auth": {site: {"ok": False, "last_check": None} for site in SITES},
        "shared_tokens": False,      # True = deliberate shared service account
    }


def read_json(path: Path) -> dict:
    """Read a JSON object; an absent, unreadable or invalid file gives ``{}``."""
    try:
        if not path.exists():
            return {}
        with path.open("r", encoding="utf-8") as handle:
            data = json.load(handle)
    except OSError as exc:
        log.warning("Config not accessible (%s: %s), ignored.", path, exc.strerror or exc)
        return {}
    except json.JSONDecodeError:
        log.warning("Unreadable config (%s), ignored.", path)
        return {}
    return data if isinstance(data, dict) else {}


def merge_config(cfg: Config, data: dict) -> None:
    """Merge ``data`` into ``cfg``, two levels deep (for key_names and auth)."""
    for key, value in data.items():
        current = cfg.get(key)
        if isinstance(value, dict) and isinstance(current, dict):
            for sub, subval in value.items():
                if isinstance(subval, dict) and isinstance(current.get(sub), dict):
                    current[sub].update(subval)
                else:
                    current[sub] = subval
        else:
            cfg[key] = value


def expand_token_dir(template: str) -> str:
    """'~/x', '{home}/x', '/data/tokens/{user}' -> absolute path for the current user."""
    text = str(template).replace("{user}", getpass.getuser()).replace("{home}", str(Path.home()))
    return str(Path(text).expanduser().resolve())


def is_personal_template(template: str) -> bool:
    """True if the template resolves to a different directory for each user."""
    text = str(template)
    return text.startswith("~") or "{user}" in text or "{home}" in text


def load_config(user_path: Path) -> Config:
    """
    Load the effective configuration.

    Only the personal part is ever written back (see :func:`save_config`); the
    shared part and its origin are kept under ``cfg["_meta"]``.
    """
    cfg = default_config()
    shared_path = shared_config_path()
    shared = read_json(shared_path)
    user = read_json(user_path)
    merge_config(cfg, {k: v for k, v in shared.items() if k in SHARED_KEYS})
    merge_config(cfg, user)
    cfg["token_dir"] = expand_token_dir(cfg.get("token_dir") or DEFAULT_TOKEN_DIR)
    cfg["_meta"] = {"shared_path": shared_path, "shared": shared, "user": user}
    return cfg


def set_user(cfg: Config, key: str, value: object) -> None:
    """Change a setting and record it as a personal one."""
    cfg[key] = value
    cfg["_meta"]["user"][key] = value


def save_config(path: Path, cfg: Config) -> None:
    """Write the personal settings only, with private permissions."""
    path.parent.mkdir(parents=True, exist_ok=True)
    restrict_permissions(path.parent, 0o700)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(cfg["_meta"]["user"], handle, indent=2, sort_keys=True)
    restrict_permissions(path, 0o600)


def restrict_permissions(path: Path, mode: int) -> None:
    """chmod ``path`` to ``mode`` if we own it and it is group/world accessible (POSIX)."""
    if os.name == "nt":
        return
    try:
        st = os.stat(str(path))
        if st.st_uid == os.getuid() and st.st_mode & 0o077:
            os.chmod(str(path), mode)
    except OSError as exc:
        log.debug("Cannot restrict permissions of %s: %s", path, exc)


def secure_token_dir(cfg: Config) -> Path:
    """
    Make sure the token directory is personal before bigsdb_downloader uses it.

    The directory holds the client secret and OAuth tokens of one BIGSdb account:
    it is created as 700 with 600 files, a directory owned by another user is
    refused, and a directory inside the (shared) script folder triggers a warning.
    Skipped when the admin explicitly enabled ``shared_tokens``.

    Raises
    ------
    BigsdbError
        If the directory cannot be created/read or belongs to another user.
    """
    tok = Path(str(cfg["token_dir"]))
    try:
        tok.mkdir(parents=True, exist_ok=True)
        owner = os.stat(str(tok)).st_uid if os.name != "nt" else None
    except OSError as exc:
        raise BigsdbError(
            "Token directory not accessible: {} ({}).\nEach user must use a personal "
            "directory, e.g.: {} setup --token-dir {}".format(
                tok, exc.strerror or exc, PROG, DEFAULT_TOKEN_DIR))
    if cfg.get("shared_tokens"):
        return tok
    if os.name != "nt":
        if owner != os.getuid():
            raise BigsdbError(
                "The token directory {} belongs to another user.\n"
                "Each user must have their own tokens: set a personal directory "
                "({} setup --token-dir {}).".format(tok, PROG, DEFAULT_TOKEN_DIR))
        restrict_permissions(tok, 0o700)
        for item in tok.rglob("*"):
            restrict_permissions(item, 0o700 if item.is_dir() else 0o600)
    script_dir = Path(__file__).resolve().parent
    if tok == script_dir or script_dir in tok.parents:
        log.warning("The token directory %s is inside the script's shared directory: "
                    "other users might be able to access it.", tok)
    return tok


def now_str() -> str:
    return datetime.datetime.now().strftime("%Y-%m-%dT%H:%M")


# --------------------------------------------------------------------------
# Backends: bigsdb_downloader (authenticated) or requests (public, pre-2025)
# --------------------------------------------------------------------------
DownloaderSpec = Union[str, Tuple[Optional[str], Optional[str]], None]


def downloader_spec(cfg: Config) -> Tuple[Optional[str], Optional[str]]:
    """(downloader path, interpreter for a .py downloader) from the config."""
    return cfg.get("downloader"), cfg.get("downloader_python")


def resolve_downloader_cmd(spec: DownloaderSpec) -> List[str]:
    """
    Turn the configured downloader into the start of a command line.

    - ``.py`` file      -> ``[python, file.py]`` (python = ``downloader_python``
      if set, otherwise the interpreter running BIGSip);
    - executable / name -> ``[executable]`` (venv, pipx, PATH).
    """
    python = None
    if isinstance(spec, (tuple, list)):
        spec, python = spec
    if not spec:
        raise BigsdbError("bigsdb_downloader path not configured. Run: {} setup".format(PROG))
    path = Path(spec).expanduser()
    if path.suffix.lower() == ".py":
        if not path.is_file():
            raise BigsdbError("bigsdb_downloader not found: {}".format(path))
        return [python or sys.executable, str(path)]
    if path.is_file():
        return [str(path)]
    found = shutil.which(str(spec))
    if found:
        return [found]
    raise BigsdbError("bigsdb_downloader not found: {}".format(spec))


def utf8_env() -> Dict[str, str]:
    """Environment forcing UTF-8 I/O in bigsdb_downloader (Windows defaults to cp1252)."""
    env = dict(os.environ)
    env.setdefault("PYTHONUTF8", "1")
    env.setdefault("PYTHONIOENCODING", "utf-8")
    return env


def format_cmd(cmd: Sequence[str]) -> str:
    """Command line as the user would type it (shell quoting of the platform)."""
    if os.name == "nt":
        return subprocess.list2cmdline(list(cmd))
    return " ".join(shlex.quote(part) for part in cmd)


def check_format(path: Path, fmt: str) -> None:
    """
    Check that a downloaded file looks like the expected format.

    An HTTP error page or a JSON error saved instead of the data would otherwise
    be kept as a valid result and silently break downstream analyses.

    Raises
    ------
    BigsdbError
        If the file is empty or does not start like a FASTA / TSV file.
    """
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        head = ""
        for line in handle:
            if line.strip():
                head = line.strip()
                break
    if not head:
        raise BigsdbError("Empty response")
    if fmt == "fasta" and not head.startswith(">"):
        raise BigsdbError("Not a FASTA file (starts with: {})".format(short_message(head, 80)))
    if fmt == "tsv" and "\t" not in head:
        raise BigsdbError("Not a TSV table (starts with: {})".format(short_message(head, 80)))


def finalize_download(tmp: Path, dest: Path, fmt: str) -> int:
    """Validate the temporary file then move it into place (atomic, no partial files)."""
    try:
        check_format(tmp, fmt)
    except BigsdbError:
        tmp.unlink()
        raise
    os.replace(str(tmp), str(dest))
    return dest.stat().st_size


class DownloaderBackend:
    """Run API requests through bigsdb_downloader (OAuth 1.0a authentication)."""

    authenticated = True

    def __init__(self, downloader: DownloaderSpec, token_dir: Optional[str], key_name: str,
                 site: str, timeout: int = DEFAULT_TIMEOUT, dry_run: bool = False) -> None:
        self.downloader = downloader
        self.token_dir = str(Path(token_dir).expanduser()) if token_dir else None
        self.key_name = key_name
        self.site = site
        self.timeout = timeout
        self.dry_run = dry_run

    def build_cmd(self, url: str, output_file: Optional[Path] = None) -> List[str]:
        cmd = resolve_downloader_cmd(self.downloader) + [
            "--key_name", self.key_name, "--site", self.site, "--url", str(url)]
        if self.token_dir:
            cmd += ["--token_dir", self.token_dir]
        # --cron: fail with a message instead of prompting for credentials, which
        # would hang forever since stdout/stderr are captured.
        cmd.append("--cron")
        if output_file is not None:
            cmd += ["--output_file", str(output_file)]
        return cmd

    def _run(self, url: str, output_file: Optional[Path] = None) -> str:
        cmd = self.build_cmd(url, output_file)
        if self.dry_run:
            print("  [dry-run] " + format_cmd(cmd))
            return ""
        log.debug("Running: %s", format_cmd(cmd))
        try:
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 universal_newlines=True, encoding="utf-8", errors="replace",
                                 timeout=self.timeout, env=utf8_env())
        except subprocess.TimeoutExpired:
            raise BigsdbError("Timeout ({}s) on {}".format(self.timeout, url))
        except OSError as exc:
            raise BigsdbError("Cannot start bigsdb_downloader: {}".format(exc))
        if res.returncode != 0:
            msg = (res.stderr or res.stdout or "").strip()
            fatal = explain_downloader_error(msg, self.site)
            if fatal:
                raise fatal
            raise BigsdbError(msg or "bigsdb_downloader failed (code {})".format(res.returncode))
        return res.stdout

    def get_text(self, url: str) -> str:
        return self._run(url)

    def save(self, url: str, dest: Path, fmt: str = "fasta") -> int:
        """Download ``url`` to ``dest`` through a ``.part`` file; return the size in bytes."""
        dest = Path(dest)
        if self.dry_run:
            self._run(url, output_file=dest)
            return 0
        dest.parent.mkdir(parents=True, exist_ok=True)
        tmp = dest.with_name(dest.name + ".part")
        if tmp.exists():
            tmp.unlink()
        self._run(url, output_file=tmp)
        if not tmp.exists():
            raise BigsdbError("Empty response")
        return finalize_download(tmp, dest, fmt)

    def run_setup(self, database: str) -> None:
        """Interactive ``bigsdb_downloader --setup`` (terminal stays attached)."""
        cmd = resolve_downloader_cmd(self.downloader) + [
            "--key_name", self.key_name, "--site", self.site, "--db", database, "--setup"]
        if self.token_dir:
            cmd += ["--token_dir", self.token_dir]
        print("\n# Setting up BIGSdb access ({}) #".format(self.site))
        print("Prerequisites: an account on {} and an API key (client id + client secret)."
              .format(SITE_WEB[self.site]))
        print("  - Pasteur: https://bigsdb.pasteur.fr/requesting-api-key/")
        print("  - PubMLST: create a personal key at https://pubmlst.org/bigsdb")
        print("Steps:")
        print("  1. enter the client id (24 characters) and the client secret (42 characters)")
        print("     -> only asked the first time for this key;")
        print("  2. open the URL shown in a web browser, log in and authorise access;")
        print("  3. paste the verification code you get.")
        print("The token is valid for every database of the site your account is registered for.")
        print("Command: " + format_cmd(cmd) + "\n")
        res = subprocess.run(cmd, env=utf8_env())
        if res.returncode != 0:
            raise BigsdbError("BIGSdb setup failed (code {}). Check the client id/secret and "
                              "the verification code, then run 'auth' again."
                              .format(res.returncode))


class PublicBackend:
    """Plain HTTP requests without authentication (data produced before 2025 only)."""

    authenticated = False

    def __init__(self, timeout: int = DEFAULT_TIMEOUT, dry_run: bool = False) -> None:
        if requests is None:
            raise BigsdbError("Public mode needs 'requests': pip install requests")
        self.timeout = timeout
        self.dry_run = dry_run
        self.session = requests.Session()
        self.session.headers["User-Agent"] = "BIGSip/{}".format(__version__)

    def get_text(self, url: str) -> str:
        if self.dry_run:
            print("  [dry-run] GET " + url)
            return ""
        try:
            response = self.session.get(url, timeout=self.timeout)
        except requests.RequestException as exc:
            raise BigsdbError(str(exc))
        if response.status_code != 200:
            raise BigsdbError("HTTP {}: {}".format(response.status_code,
                                                   short_message(response.text)))
        return response.text

    def save(self, url: str, dest: Path, fmt: str = "fasta") -> int:
        dest = Path(dest)
        text = self.get_text(url)
        if self.dry_run:
            return 0
        dest.parent.mkdir(parents=True, exist_ok=True)
        tmp = dest.with_name(dest.name + ".part")
        tmp.write_text(text, encoding="utf-8")
        return finalize_download(tmp, dest, fmt)


Backend = Union[DownloaderBackend, PublicBackend]


# --------------------------------------------------------------------------
# BIGSdb client: catalogue, schemes and data URLs
# --------------------------------------------------------------------------
def organism_of(name: str) -> str:
    """'pubmlst_klebsiella_seqdef' -> 'klebsiella'; 'Klebsiella' -> 'klebsiella'."""
    name = name.strip().lower().replace(" ", "_").replace("-", "_")
    for suffix in ("_" + DbKind.SEQDEF.value, "_" + DbKind.ISOLATES.value):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    if name.startswith("pubmlst_"):
        name = name[len("pubmlst_"):]
    return name


class BigsdbFetcher:
    """
    Read the site catalogue and build the data URLs.

    Metadata (databases, schemes, loci) are public and read over plain HTTP when
    ``requests`` is available; data downloads always go through ``backend``.
    """

    def __init__(self, site: str, backend: Backend) -> None:
        self.site = site
        self.api = SITES[site]
        self.backend = backend
        self._public: Optional[PublicBackend] = None
        self._organisms: Optional[Dict[str, Organism]] = None

    # ---- metadata -----------------------------------------------------------
    def _meta_json(self, url: str):
        if requests is not None:
            try:
                if self._public is None:
                    self._public = PublicBackend(timeout=METADATA_TIMEOUT)
                return json.loads(self._public.get_text(url))
            except Exception as exc:  # network, TLS, proxy or JSON error
                log.debug("Public metadata read failed (%s), using bigsdb_downloader.", exc)
        # Metadata must really be read, even in dry-run mode.
        dry_run = getattr(self.backend, "dry_run", False)
        self.backend.dry_run = False
        try:
            return json.loads(self.backend.get_text(url))
        finally:
            self.backend.dry_run = dry_run

    def list_organisms(self) -> Dict[str, Organism]:
        """
        Group the site databases by organism.

        Route: ``GET /db``. Database names follow ``pubmlst_<organism>_<kind>``;
        the organism key is the name stripped of the prefix and the kind suffix.
        """
        if self._organisms is None:
            organisms: Dict[str, Organism] = {}
            for group in self._meta_json("{}/db".format(self.api)):
                for db in group.get("databases", []):
                    name = db.get("name") or ""
                    kind = next((k for k in DbKind if name.endswith("_" + k.value)), None)
                    if kind is None:
                        continue
                    org = organisms.setdefault(organism_of(name), Organism(organism_of(name)))
                    setattr(org, kind.value, name)
                    if not org.description:
                        org.description = (group.get("description") or group.get("name")
                                           or db.get("description", ""))
            self._organisms = dict(sorted(organisms.items()))
        return self._organisms

    def list_schemes(self, database: str, with_pk: bool = False) -> List[Scheme]:
        """
        Schemes of a seqdef database, sorted by ID.

        Route: ``GET /db/{database}/schemes``; ``?with_pk=1`` restricts the list to
        schemes with a primary key field, i.e. the only ones with profiles.
        """
        url = "{}/db/{}/schemes".format(self.api, database) + ("?with_pk=1" if with_pk else "")
        schemes = [Scheme(entry["scheme"].rstrip("/").split("/")[-1],
                          entry.get("description", ""), True if with_pk else None)
                   for entry in self._meta_json(url).get("schemes", [])]
        return sorted(schemes, key=lambda s: int(s.id) if s.id.isdigit() else 0)

    def list_schemes_with_profiles(self, database: str) -> List[Scheme]:
        """All schemes, with ``has_profiles`` set from the ``with_pk`` list when available."""
        schemes = self.list_schemes(database)
        try:
            pk_ids = {s.id for s in self.list_schemes(database, with_pk=True)}
        except Exception as exc:
            log.debug("Primary-key scheme list unavailable: %s", exc)
            return schemes
        for scheme in schemes:
            scheme.has_profiles = scheme.id in pk_ids
        return schemes

    def scheme_info(self, database: str, scheme_id: str) -> SchemeInfo:
        """Route: ``GET /db/{database}/schemes/{scheme_id}`` (loci are given as URLs)."""
        data = self._meta_json("{}/db/{}/schemes/{}".format(self.api, database, scheme_id))
        return SchemeInfo(str(scheme_id),
                          data.get("description", "scheme_{}".format(scheme_id)),
                          [u.rstrip("/").split("/")[-1] for u in data.get("loci", [])])

    # ---- data URLs ----------------------------------------------------------
    def url_genome(self, database: str, isolate_id: int, header: Optional[str] = None) -> str:
        """Route: ``GET /db/{database}/isolates/{isolate_id}/contigs_fasta[?header=]``."""
        url = "{}/db/{}/isolates/{}/contigs_fasta".format(self.api, database, isolate_id)
        return url + ("?header={}".format(header) if header else "")

    def url_locus_fasta(self, database: str, locus: str) -> str:
        """Route: ``GET /db/{database}/loci/{locus}/alleles_fasta``."""
        return "{}/db/{}/loci/{}/alleles_fasta".format(self.api, database, quote(locus, safe=""))

    def url_profiles(self, database: str, scheme_id: str) -> str:
        """Route: ``GET /db/{database}/schemes/{scheme_id}/profiles_csv`` (tab-separated)."""
        return "{}/db/{}/schemes/{}/profiles_csv".format(self.api, database, scheme_id)


# --------------------------------------------------------------------------
# User input parsing
# --------------------------------------------------------------------------
def split_items(text: Optional[str]) -> List[str]:
    """'a, b;c d' -> ['a', 'b', 'c', 'd'] (comma, semicolon or whitespace separated)."""
    return [x for x in re.split(r"[,;\s]+", text or "") if x]


def split_names(text: Optional[str]) -> List[str]:
    """Locus names: comma, semicolon or newline separated (spaces are kept)."""
    return [x.strip() for x in re.split(r"[,;\n]+", text or "") if x.strip()]


def read_list_file(path: Union[str, Path]) -> List[str]:
    """One item per line (commas tolerated); empty lines and ``#`` comments ignored."""
    path = Path(path).expanduser()
    if not path.is_file():
        raise BigsdbError("File not found: {}".format(path))
    items: List[str] = []
    with path.open("r", encoding="utf-8-sig") as handle:
        for line in handle:
            line = line.split("#", 1)[0].strip()
            if line:
                items.extend(split_names(line))
    return items


def names_or_file(value: Optional[str], numeric: bool = False) -> List[str]:
    """An option value is either a comma-separated list or the path of a list file."""
    if value and Path(value).expanduser().is_file():
        return read_list_file(value)
    return split_items(value) if numeric else split_names(value)


def expand_ids(tokens: Iterable[str]) -> List[int]:
    """
    Expand isolate IDs and ranges: ['1', '5-8', '12'] -> [1, 5, 6, 7, 8, 12].

    Duplicates are removed, order is kept.

    Raises
    ------
    BigsdbError
        If a token is not a number or a range (e.g. a strain name).
    """
    ids: List[int] = []
    for tok in tokens:
        match = re.fullmatch(r"(\d+)\s*-\s*(\d+)", tok)
        if match:
            low, high = sorted((int(match.group(1)), int(match.group(2))))
            ids.extend(range(low, high + 1))
        elif tok.isdigit():
            ids.append(int(tok))
        else:
            raise BigsdbError(
                "'{}' is not an isolate ID. Genomes are downloaded by NUMERIC ID "
                "(e.g. 1,2,10-20), not by strain name: the BIGSdb URL is built from the ID."
                .format(tok))
    return unique(ids)


def collect_scheme_ids(values: Optional[Sequence[str]]) -> List[str]:
    """Scheme IDs from repeated/comma-separated options; names are refused."""
    ids = [item for value in values or [] for item in split_items(value)]
    for item in ids:
        if not item.isdigit():
            raise BigsdbError(
                "'{}' is not a scheme ID. Use the NUMERIC scheme ID (e.g. 4), not its name; "
                "list: {} schemes --organism <org>".format(item, PROG))
    return unique(ids)


def unique(seq: Iterable) -> list:
    return list(dict.fromkeys(seq))


def safe_name(text: str, max_len: int = 80) -> str:
    """File-system safe name (scheme descriptions and locus names can contain anything)."""
    text = re.sub(r"[^\w.\-]+", "_", str(text)).strip("_")
    return text[:max_len] or "unnamed"


# --------------------------------------------------------------------------
# Organism -> database, scheme checks
# --------------------------------------------------------------------------
def resolve_database(fetcher: BigsdbFetcher, organism: str, kind: DbKind) -> str:
    """
    Build and check the database name for an organism and a data type.

    genomes -> ``<organism>_isolates``; loci / profiles -> ``<organism>_seqdef``.
    A full database name is accepted too (prefix and suffix are ignored).

    Raises
    ------
    BigsdbError
        Unknown organism (with close names), or organism without this database type.
    """
    org = organism_of(organism)
    if not org:
        raise BigsdbError("Empty organism name.")
    try:
        organisms = fetcher.list_organisms()
    except Exception as exc:
        db = "pubmlst_{}_{}".format(org, kind.value)
        log.warning("Organism list unavailable (%s); using '%s' without checking.",
                    short_message(str(exc), 80), db)
        return db

    if org not in organisms:
        close = difflib.get_close_matches(org, list(organisms), n=5, cutoff=0.5)
        close += [name for name in organisms if org in name and name not in close][:5]
        msg = "Unknown organism '{}' on {}.".format(organism, fetcher.site)
        if close:
            msg += " Did you mean: " + ", ".join(close) + "?"
        raise BigsdbError(msg + "\nList of organisms: {} organisms --site {}"
                          .format(PROG, fetcher.site))

    db = organisms[org].db(kind)
    if not db:
        available = [k.label for k in DbKind if organisms[org].db(k)]
        raise BigsdbError(
            "'{}' has no {} database on {}: cannot download {}.\n"
            "Available for this organism: {}.".format(
                org, kind.value, fetcher.site, kind.data, ", ".join(available) or "nothing"))
    return db


def schemes_table(schemes: Sequence[Scheme], indent: str = "  ") -> str:
    lines = [indent + "{:<6} | {}".format("ID", "Description"), indent + "-" * 58]
    lines += [indent + "{:<6} | {}".format(s.id, s.description) for s in schemes]
    return "\n".join(lines)


def validate_schemes(fetcher: BigsdbFetcher, database: str, org: str,
                     scheme_ids: Sequence[str]) -> None:
    """Raise BigsdbError listing the available schemes if an ID does not exist."""
    try:
        schemes = fetcher.list_schemes(database)
    except Exception as exc:
        log.debug("Scheme list unavailable, IDs not checked: %s", exc)
        return
    known = {s.id for s in schemes}
    missing = [sid for sid in scheme_ids if str(sid) not in known]
    if missing:
        raise BigsdbError(
            "Scheme(s) {} not found for {}.\nAvailable schemes:\n{}\n"
            "(also listed by: {} schemes --organism {})".format(
                ", ".join(map(str, missing)), org, schemes_table(schemes), PROG, org))


def check_profile_schemes(fetcher: BigsdbFetcher, database: str, org: str,
                          scheme_ids: Sequence[str]) -> None:
    """
    Refuse schemes without a primary key before requesting their profiles.

    A profile table maps an identifier (ST, cgST...) to an allele combination;
    schemes that are only a set of loci (no primary key) have none, and
    ``profiles_csv`` returns HTTP 400 for them.
    """
    try:
        descriptions = {s.id: s.description for s in fetcher.list_schemes(database)}
        with_pk = fetcher.list_schemes(database, with_pk=True)
    except Exception as exc:
        log.debug("Primary-key check skipped: %s", exc)
        return
    pk_ids = {s.id for s in with_pk}
    no_pk = [sid for sid in scheme_ids if str(sid) not in pk_ids]
    if no_pk:
        raise BigsdbError(
            "Scheme(s) without profiles: {}.\n"
            "A scheme only has profiles if it has a primary key (e.g. ST for an MLST); "
            "this one is just a set of loci.\n"
            "  -> For its sequences: {} loci --organism {} --scheme {}\n"
            "Schemes of {} with profiles:\n{}".format(
                ", ".join("{} ({})".format(sid, descriptions.get(str(sid), "?")) for sid in no_pk),
                PROG, org, ",".join(map(str, no_pk)), org,
                schemes_table(with_pk) if with_pk else "  (none)"))


# --------------------------------------------------------------------------
# Download engine: resume, parallelism, failure list, provenance log
# --------------------------------------------------------------------------
def human_size(n_bytes: float) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n_bytes < 1024:
            return "{:.0f} {}".format(n_bytes, unit)
        n_bytes /= 1024.0
    return "{:.1f} TB".format(n_bytes)


def run_jobs(fetcher: BigsdbFetcher, jobs: Sequence[Job], overwrite: bool = False,
             n_workers: int = 1, label: str = "item") -> BatchResult:
    """
    Download ``jobs``, skipping files already present unless ``overwrite``.

    The first download always runs alone: it validates authentication (a
    FatalError stops everything at once) and lets bigsdb_downloader renew the
    session token before several processes would rewrite its token file.
    """
    result = BatchResult()
    todo: List[Job] = []
    for job in jobs:
        if not overwrite and job.dest.is_file() and job.dest.stat().st_size > 0:
            result.skipped.append(job.key)
        else:
            todo.append(job)
    if result.skipped:
        log.info("%d %s(s) already present, skipped (--overwrite to force).",
                 len(result.skipped), label)

    dry_run = getattr(fetcher.backend, "dry_run", False)
    total = len(todo)

    def download(job: Job) -> int:
        return fetcher.backend.save(job.url, job.dest, job.fmt)

    def report(index: int, job: Job, size: Optional[int], error: Optional[str]) -> None:
        if error is not None:
            result.failed.append((job.key, error))
            status = "FAILED - " + short_message(error)
        else:
            result.ok.append(job.key)
            status = "dry-run" if dry_run else "OK ({})".format(human_size(size or 0))
        log.info("[%d/%d] %s %s: %s", index, total, label, job.key, status)

    n_sequential = total if n_workers <= 1 else (1 if fetcher.backend.authenticated else 0)
    for index, job in enumerate(todo[:n_sequential], 1):
        try:
            report(index, job, download(job), None)
        except FatalError:
            raise  # the whole batch would fail the same way
        except BigsdbError as exc:
            report(index, job, None, str(exc))

    remaining = todo[n_sequential:]
    if remaining:
        with ThreadPoolExecutor(max_workers=n_workers) as pool:
            futures = {pool.submit(download, job): job for job in remaining}
            for index, future in enumerate(as_completed(futures), n_sequential + 1):
                job = futures[future]
                try:
                    report(index, job, future.result(), None)
                except FatalError:
                    for pending in futures:
                        pending.cancel()
                    raise
                except BigsdbError as exc:
                    report(index, job, None, str(exc))
    return result


def write_failures(outdir: Path, kind: str, failed: Sequence[Tuple[str, str]]) -> Optional[Path]:
    """Write the failed items, one per line, directly reusable as --ids / --loci input."""
    if not failed:
        return None
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    path = outdir / "failed_{}_{}.txt".format(kind, stamp)
    path.parent.mkdir(parents=True, exist_ok=True)
    option = {"genomes": "--ids", "loci": "--loci"}.get(kind, "--scheme")
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# Failed {} - retry with: {} {}\n".format(kind, option, path.name))
        for key, msg in failed:
            handle.write("{}\t# {}\n".format(key, short_message(msg, 150)))
    return path


def write_download_log(outdir: Path, fetcher: BigsdbFetcher, database: str, data_type: str,
                       request: str, n_requested: int, result: BatchResult) -> None:
    """
    Append one line per run to ``download_log.tsv`` for data provenance.

    Records when, from where and with which request the files of the folder
    were obtained, so an analysis can cite the exact data source later.
    """
    path = outdir / DOWNLOAD_LOG
    path.parent.mkdir(parents=True, exist_ok=True)
    new_file = not path.exists()
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        if new_file:
            writer.writerow(["date", "bigsip_version", "site", "database", "data_type",
                             "request", "authenticated", "requested", "downloaded",
                             "already_present", "failed"])
        writer.writerow([datetime.datetime.now().isoformat(timespec="seconds"), __version__,
                         fetcher.site, database, data_type, request,
                         "yes" if fetcher.backend.authenticated else "no (public, pre-2025)",
                         n_requested, len(result.ok), len(result.skipped), len(result.failed)])


def summarize(kind: str, result: BatchResult, outdir: Path, dry_run: bool = False) -> None:
    """Print the run summary and write the failure list."""
    if dry_run:
        print("\n[dry-run] No request was run.")
        return
    print("\n=== {} summary ===".format(kind.capitalize()))
    print("  downloaded      : {}".format(len(result.ok)))
    print("  already present : {}".format(len(result.skipped)))
    print("  failed          : {}".format(len(result.failed)))
    print("  directory       : {}".format(outdir.resolve()))
    failures = write_failures(outdir, kind, result.failed)
    if failures:
        print("  failure list    : {}".format(failures))


def describe_ids(ids: Sequence[int]) -> str:
    """Compact form of an ID list for the provenance log: [1, 2, 3, 7] -> '1-3,7'."""
    runs: List[List[int]] = []
    for value in sorted(set(ids)):
        if runs and value == runs[-1][1] + 1:
            runs[-1][1] = value
        else:
            runs.append([value, value])
    return ",".join(str(a) if a == b else "{}-{}".format(a, b) for a, b in runs)


# --------------------------------------------------------------------------
# Download actions
# --------------------------------------------------------------------------
def do_genomes(fetcher: BigsdbFetcher, organism: str, ids: Sequence[int], outdir: str,
               header: Optional[str] = None, merge: bool = False, overwrite: bool = False,
               jobs: int = 1, dry_run: bool = False) -> bool:
    """Download the contigs of each isolate into ``<outdir>/<organism>/genomes/<id>.fasta``."""
    database = resolve_database(fetcher, organism, DbKind.ISOLATES)
    org = organism_of(organism)
    org_dir = Path(outdir) / org
    target = org_dir / "genomes"
    log.info("%d genome(s) %s (database %s) -> %s", len(ids), org, database, target)
    job_list = [Job(str(i), fetcher.url_genome(database, i, header),
                    target / "{}.fasta".format(i)) for i in ids]
    result = run_jobs(fetcher, job_list, overwrite, jobs, "isolate")
    if merge and not dry_run:
        merged = org_dir / "{}_genomes.fasta".format(org)
        with merged.open("w", encoding="utf-8") as out:
            for i in ids:
                fasta = target / "{}.fasta".format(i)
                if fasta.is_file():
                    out.write(fasta.read_text(encoding="utf-8").rstrip("\n") + "\n")
        log.info("Merged multi-FASTA: %s", merged)
    if not dry_run:
        write_download_log(org_dir, fetcher, database, "genomes",
                           "ids=" + describe_ids(ids), len(job_list), result)
    summarize("genomes", result, org_dir, dry_run)
    return not result.failed


def do_loci(fetcher: BigsdbFetcher, organism: str, loci: Optional[Sequence[str]] = None,
            scheme_ids: Optional[Sequence[str]] = None, outdir: str = DEFAULT_OUTDIR,
            overwrite: bool = False, jobs: int = 1, dry_run: bool = False) -> bool:
    """
    Download allele FASTA files.

    Chosen loci go to ``<organism>/alleles/``; every locus of a scheme goes to
    ``<organism>/scheme_<id>_<description>/alleles/`` with the locus list in
    ``loci.txt``.
    """
    database = resolve_database(fetcher, organism, DbKind.SEQDEF)
    org = organism_of(organism)
    org_dir = Path(outdir) / org
    if scheme_ids:
        validate_schemes(fetcher, database, org, scheme_ids)

    groups: List[Tuple[Path, List[str], str]] = []   # (folder, loci, provenance request)
    if loci:
        groups.append((org_dir / "alleles", unique(loci), "loci=" + ",".join(unique(loci))))
    for sid in scheme_ids or []:
        info = fetcher.scheme_info(database, sid)
        folder = org_dir / "scheme_{}_{}".format(sid, safe_name(info.description))
        log.info("Scheme %s (%s): %d loci", sid, info.description, len(info.loci))
        if not info.loci:
            log.warning("No locus found for scheme %s.", sid)
            continue
        groups.append((folder / "alleles", info.loci, "scheme={}".format(sid)))
        if not dry_run:
            folder.mkdir(parents=True, exist_ok=True)
            (folder / "loci.txt").write_text("\n".join(info.loci) + "\n", encoding="utf-8")

    all_ok = True
    for folder, locus_list, request in groups:
        log.info("%d loci %s (database %s) -> %s", len(locus_list), org, database, folder)
        job_list = [Job(locus, fetcher.url_locus_fasta(database, locus),
                        folder / "{}.fasta".format(safe_name(locus, 150))) for locus in locus_list]
        result = run_jobs(fetcher, job_list, overwrite, jobs, "locus")
        if not dry_run:
            write_download_log(org_dir, fetcher, database, "loci", request, len(job_list), result)
        summarize("loci", result, folder.parent, dry_run)
        all_ok = all_ok and not result.failed
    return all_ok


def do_profiles(fetcher: BigsdbFetcher, organism: str, scheme_ids: Sequence[str], outdir: str,
                overwrite: bool = False, dry_run: bool = False) -> bool:
    """Download ``profiles_csv`` of each scheme into ``scheme_<id>_<name>/profiles.tsv``."""
    database = resolve_database(fetcher, organism, DbKind.SEQDEF)
    org = organism_of(organism)
    org_dir = Path(outdir) / org
    validate_schemes(fetcher, database, org, scheme_ids)
    check_profile_schemes(fetcher, database, org, scheme_ids)
    job_list = []
    for sid in scheme_ids:
        try:
            description = fetcher.scheme_info(database, sid).description
        except Exception as exc:
            log.debug("Description of scheme %s unavailable: %s", sid, exc)
            description = "scheme_{}".format(sid)
        folder = org_dir / "scheme_{}_{}".format(sid, safe_name(description))
        job_list.append(Job(str(sid), fetcher.url_profiles(database, sid),
                            folder / "profiles.tsv", fmt="tsv"))
    log.info("%d scheme(s): profiles %s (database %s)", len(job_list), org, database)
    result = run_jobs(fetcher, job_list, overwrite, 1, "scheme")
    if not dry_run:
        write_download_log(org_dir, fetcher, database, "profiles",
                           "schemes=" + ",".join(map(str, scheme_ids)), len(job_list), result)
    summarize("profiles", result, org_dir, dry_run)
    return not result.failed


# --------------------------------------------------------------------------
# Display
# --------------------------------------------------------------------------
def filter_organisms(organisms: Dict[str, Organism], search: Optional[str] = None,
                     need: Optional[DbKind] = None) -> List[Organism]:
    """Organisms matching ``search`` (name or description) and having ``need``."""
    rows = []
    for org in organisms.values():
        if need and not org.db(need):
            continue
        if search and search.lower() not in (org.name + " " + org.description).lower():
            continue
        rows.append(org)
    return rows


def print_organisms(organisms: Dict[str, Organism], search: Optional[str] = None,
                    need: Optional[DbKind] = None, site: str = "") -> List[Organism]:
    rows = filter_organisms(organisms, search, need)
    if not rows:
        print("No matching organism.")
        return rows
    width = max(12, max(len(org.name) for org in rows))
    print("\nAvailable organisms{}:\n".format(" on " + site if site else ""))
    print("{:<{w}}  {:<7}  {:<13}  {}".format("Organism", "Genomes", "Loci/profiles",
                                             "Description", w=width))
    print("-" * (width + 60))
    for org in rows:
        print("{:<{w}}  {:<7}  {:<13}  {}".format(
            org.name, "yes" if org.isolates else "-", "yes" if org.seqdef else "-",
            org.description, w=width))
    print("\n{} organism(s). Use the name in the first column with --organism.".format(len(rows)))
    return rows


def print_schemes(schemes: Sequence[Scheme], organism: str) -> None:
    if not schemes:
        print("No scheme for {}.".format(organism))
        return
    print("\nSchemes of {}:\n".format(organism))
    if any(s.has_profiles is not None for s in schemes):
        print("{:<6} | {:<8} | {}".format("ID", "Profiles", "Description"))
        print("-" * 60)
        for s in schemes:
            print("{:<6} | {:<8} | {}".format(s.id, "yes" if s.has_profiles else "-",
                                               s.description))
        if not all(s.has_profiles for s in schemes):
            print("\n'-' = scheme without a primary key (no ST): no profiles, "
                  "but its loci can still be downloaded with 'loci'.")
    else:
        print(schemes_table(schemes, indent=""))


def print_profile_schemes(fetcher: BigsdbFetcher, database: str,
                          org: str) -> Optional[List[Scheme]]:
    """
    Show the schemes usable with 'profiles' and name the excluded ones.

    Returns the schemes with profiles, or None if this cannot be determined.
    """
    try:
        all_schemes = fetcher.list_schemes(database)
        with_pk = fetcher.list_schemes(database, with_pk=True)
    except Exception as exc:
        log.debug("Primary-key scheme list unavailable: %s", exc)
        return None
    pk_ids = {s.id for s in with_pk}
    excluded = [s for s in all_schemes if s.id not in pk_ids]
    print("\nSchemes of {} WITH profiles (primary key, e.g. ST / cgST):\n".format(org))
    print(schemes_table(with_pk, indent="") if with_pk else "  none.")
    if excluded:
        print("\nWithout profiles (no primary key) - only through 'loci --scheme': {}".format(
            ", ".join("{} {}".format(s.id, s.description) for s in excluded)))
    return with_pk


# --------------------------------------------------------------------------
# Building fetchers from the config
# --------------------------------------------------------------------------
def print_auth_help(site: str) -> None:
    print("\n# BIGSdb authentication not configured for {site} #\n"
          "Data produced after 01/01/2025 require authenticated access to the\n"
          "BIGSdb API through bigsdb_downloader (personal API key).\n"
          "  - Pasteur: https://bigsdb.pasteur.fr/requesting-api-key/\n"
          "  - PubMLST: https://pubmlst.org/bigsdb (personal key)\n"
          "Setup:\n"
          "  {prog} setup\n"
          "  {prog} auth --site {site} --organism <organism>\n"
          "Without authentication, use --public (pre-2025 data only).\n"
          .format(site=site, prog=PROG))


def downloader_backend(cfg: Config, site: str, timeout: int,
                       dry_run: bool = False) -> DownloaderBackend:
    return DownloaderBackend(downloader_spec(cfg), cfg.get("token_dir"),
                             cfg["key_names"].get(site, site), site, timeout, dry_run)


def make_fetcher(args: argparse.Namespace, cfg: Config) -> BigsdbFetcher:
    """Fetcher used for downloads: bigsdb_downloader, or plain HTTP with --public."""
    site = args.site or cfg.get("default_site") or DEFAULT_SITE
    if args.public:
        log.info("Public mode: only data produced before 01/01/2025 are available.")
        return BigsdbFetcher(site, PublicBackend(args.timeout, args.dry_run))
    if not cfg.get("downloader"):
        print_auth_help(site)
        raise BigsdbError("bigsdb_downloader not configured.")
    secure_token_dir(cfg)
    if not cfg["auth"].get(site, {}).get("ok"):
        log.warning("%s authentication has not been validated for %s ('auth' command). "
                    "Trying anyway...", site, getpass.getuser())
    return BigsdbFetcher(site, downloader_backend(cfg, site, args.timeout, args.dry_run))


def meta_fetcher(args: argparse.Namespace, cfg: Config) -> BigsdbFetcher:
    """Fetcher used for listings only: no authentication needed, a missing config is fine."""
    site = args.site or cfg.get("default_site") or DEFAULT_SITE
    if cfg.get("downloader") and not args.public:
        return BigsdbFetcher(site, downloader_backend(cfg, site, args.timeout))
    return BigsdbFetcher(site, PublicBackend(args.timeout))


# --------------------------------------------------------------------------
# Interactive helpers
# --------------------------------------------------------------------------
def ask(prompt: str, default: Optional[str] = None) -> str:
    """Prompt the user; on end of input, return ``default`` or stop if there is none."""
    suffix = " [{}]".format(default) if default not in (None, "") else ""
    try:
        value = input("{}{}: ".format(prompt, suffix)).strip()
    except EOFError:
        if default is None:
            raise BigsdbError("Input interrupted.")
        value = ""
    return value or (default if default is not None else "")


def ask_yes(prompt: str, default: bool = True) -> bool:
    answer = ask("{} ({})".format(prompt, "Y/n" if default else "y/N"), "").lower()
    return default if not answer else answer in ("y", "yes", "o", "oui")


def ask_until(prompt: str, parse: Callable[[str], object]):
    """Ask again while the input is empty or ``parse`` raises BigsdbError."""
    while True:
        value = ask(prompt, "")
        if not value:
            print("[ERROR] Input required.")
            continue
        try:
            result = parse(value)
            if result:
                return result
            print("[ERROR] No value recognised.")
        except BigsdbError as exc:
            print("[ERROR] {}".format(exc))


def choose_from(rows: Sequence, label: Callable[[object], str], prompt: str,
                key: Callable[[object], str] = lambda row: getattr(row, "name", "")):
    """Numbered menu; the exact name of an entry is accepted as well as its number."""
    for index, row in enumerate(rows, 1):
        print("  {:>3}. {}".format(index, label(row)))
    while True:
        value = ask(prompt)
        if value.isdigit() and 1 <= int(value) <= len(rows):
            return rows[int(value) - 1]
        for row in rows:
            if value and value == key(row):
                return row
        print("Invalid choice.")


def pick_organism(fetcher: BigsdbFetcher, kind: DbKind) -> str:
    """
    Ask for ONE organism having a ``kind`` database.

    The input can be the exact name, or part of a name/description used to
    filter the list. An empty input shows the list and asks again: there is
    deliberately no 'all organisms' choice before a download.
    """
    organisms = fetcher.list_organisms()
    while True:
        value = ask("Organism (name, or part of the name to filter the list)", "").strip()
        if not value:
            print("\n[ERROR] Missing organism.")
            print_organisms(organisms, need=kind, site=fetcher.site)
            continue
        org = organism_of(value)
        if org in organisms and organisms[org].db(kind):
            return org

        rows = filter_organisms(organisms, value, kind)
        if not rows:
            other = filter_organisms(organisms, value)
            if other:
                print("\n[ERROR] '{}' matches: {}; but none of them has {} on {}.".format(
                    value, ", ".join(o.name for o in other), kind.label, fetcher.site))
                continue
            try:
                resolve_database(fetcher, value, kind)
            except BigsdbError as exc:
                print("\n[ERROR] {}".format(exc))
            continue
        if len(rows) == 1:
            if ask_yes("-> {} ({})?".format(rows[0].name, rows[0].description)):
                return rows[0].name
            continue
        print("\n{} organisms match '{}':".format(len(rows), value))
        return choose_from(rows, lambda o: "{:<20} {}".format(o.name, o.description),
                           "Number or name of the organism").name


# --------------------------------------------------------------------------
# Commands: setup, auth, status
# --------------------------------------------------------------------------
def find_downloader() -> Optional[str]:
    """Suggest 'bigsdb-downloader' next to the running Python (same venv), else on the PATH."""
    bindir = Path(sys.executable).parent
    for name in ("bigsdb-downloader", "bigsdb-downloader.exe"):
        if (bindir / name).is_file():
            return str(bindir / name)
    return shutil.which("bigsdb-downloader")


def check_downloader(downloader: str, python: Optional[str] = None) -> str:
    """Run ``--help`` on the downloader, report problems, return its absolute path."""
    try:
        cmd = resolve_downloader_cmd((downloader, python))
    except BigsdbError as exc:
        raise BigsdbError(str(exc) + "\nInstall: pip install bigsdb-downloader  "
                          "(https://github.com/kjolley/BIGSdb_downloader)")
    res = subprocess.run(cmd + ["--help"], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         universal_newlines=True, encoding="utf-8", errors="replace",
                         env=utf8_env())
    if res.returncode == 0 and "--key_name" in res.stdout:
        print("[OK] bigsdb_downloader works ({}).".format(format_cmd(cmd)))
    else:
        broken = explain_downloader_error(res.stderr or res.stdout, DEFAULT_SITE)
        log.warning("'%s --help' failed.", format_cmd(cmd))
        print(str(broken) if broken else (res.stderr or res.stdout)[:400])
    path = Path(downloader).expanduser()
    return str(path.resolve()) if path.exists() else downloader


def _setup_downloader(args: argparse.Namespace, cfg: Config, interactive: bool) -> None:
    """Downloader path: CLI option > shared config > existing personal setting > prompt."""
    shared = cfg["_meta"]["shared"]
    if args.python:
        set_user(cfg, "downloader_python", args.python)
    python = cfg.get("downloader_python")
    if args.downloader:
        set_user(cfg, "downloader", check_downloader(args.downloader, python))
    elif shared.get("downloader"):
        print("\nbigsdb_downloader (shared): {}".format(shared["downloader"]))
        check_downloader(shared["downloader"], python)
    elif cfg.get("downloader") and not interactive:
        check_downloader(cfg["downloader"], python)
    else:
        downloader = None
        while not downloader:
            print("\nPath to bigsdb_downloader: either the bigsdb_downloader.py script,\n"
                  "or the 'bigsdb-downloader' executable (pip/pipx install bigsdb-downloader).")
            downloader = ask("Path", cfg.get("downloader") or find_downloader())
        set_user(cfg, "downloader", check_downloader(downloader, python))


def _setup_token_dir(args: argparse.Namespace, cfg: Config, interactive: bool) -> None:
    token_dir = args.token_dir
    if not token_dir and interactive:
        print("\nPERSONAL directory where your BIGSdb keys and tokens are stored.")
        token_dir = ask("Token directory", cfg["token_dir"])
    if token_dir:
        cfg["_meta"]["user"]["token_dir"] = token_dir
        cfg["token_dir"] = expand_token_dir(token_dir)
    secure_token_dir(cfg)


def _setup_site_and_key(args: argparse.Namespace, cfg: Config, interactive: bool) -> str:
    current = cfg.get("default_site", DEFAULT_SITE)
    site = args.site or current
    if interactive and not args.site:
        site = ask("Default site (Pasteur/PubMLST)", current)
    if site not in SITES:
        raise BigsdbError("Unknown site: {} (choices: {})".format(site, ", ".join(SITES)))
    if site != current or args.site:
        set_user(cfg, "default_site", site)
    key = args.key_name or (ask("Name of your API key for {}".format(site),
                                cfg["key_names"].get(site, site)) if interactive else None)
    if key:
        cfg["key_names"][site] = key
        set_user(cfg, "key_names", cfg["key_names"])
    return site


def cmd_setup(args: argparse.Namespace, cfg: Config, cfg_file: Path) -> int:
    """Personal setup (or shared setup with --shared)."""
    if args.shared:
        return cmd_setup_shared(args, cfg)
    print("\n# Personal BIGSip configuration ({}) #".format(getpass.getuser()))
    print("File: {}".format(cfg_file))
    if cfg["_meta"]["shared"]:
        print("Shared config found: {}".format(cfg["_meta"]["shared_path"]))
    interactive = not (args.downloader or args.token_dir or args.key_name or args.site
                       or args.python)
    _setup_downloader(args, cfg, interactive)
    _setup_token_dir(args, cfg, interactive)
    site = _setup_site_and_key(args, cfg, interactive)
    save_config(cfg_file, cfg)
    print("\n[OK] Personal configuration saved in {}".format(cfg_file))
    print("     Tokens: {}".format(cfg["token_dir"]))
    print("Next step: {} auth --site {} --organism <organism>".format(PROG, site))
    return 0


def cmd_setup_shared(args: argparse.Namespace, cfg: Config) -> int:
    """
    Shared install (admin): write ``bigsip_shared.json`` next to the script.

    It holds the downloader path and the token directory TEMPLATE, never a token.
    """
    path = cfg["_meta"]["shared_path"]
    shared = dict(cfg["_meta"]["shared"])
    print("\n# SHARED BIGSip configuration #")
    print("File: {} (read by every user)".format(path))

    downloader = args.downloader or ask("Path to bigsdb_downloader (readable by everyone)",
                                        shared.get("downloader") or find_downloader())
    if not downloader:
        raise BigsdbError("Path to bigsdb_downloader required.")
    if args.python:
        shared["downloader_python"] = args.python
    shared["downloader"] = check_downloader(downloader, shared.get("downloader_python"))

    template = args.token_dir or ask(
        "Token directory template ({user}/{home}/~ = different for each user)",
        shared.get("token_dir", DEFAULT_TOKEN_DIR))
    if not is_personal_template(template) and not args.shared_tokens:
        raise BigsdbError(
            "'{}' would be the SAME directory for everyone: users would use each other's "
            "tokens.\nUse ~, {{home}} or {{user}} (e.g. /data/bigsdb_tokens/{{user}}), "
            "or --shared-tokens for a deliberate shared service account.".format(template))
    shared["token_dir"] = template
    shared["shared_tokens"] = bool(args.shared_tokens)
    if args.site:
        shared["default_site"] = args.site
    if args.key_name:
        shared.setdefault("key_names", {})[args.site or DEFAULT_SITE] = args.key_name

    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8") as handle:
            json.dump(shared, handle, indent=2, sort_keys=True)
        if os.name != "nt":
            os.chmod(str(path), 0o644)  # readable by everyone, writable by the admin only
    except OSError as exc:
        raise BigsdbError("Cannot write {}: {}".format(path, exc))
    print("\n[OK] Shared config saved: {}".format(path))
    print("Each user now only needs to run:")
    print("  python {} auth --organism <organism>".format(Path(__file__).resolve()))
    return 0


def cmd_auth(args: argparse.Namespace, cfg: Config, cfg_file: Path) -> int:
    """Run the interactive OAuth setup of bigsdb_downloader, then test the access."""
    site = args.site or cfg.get("default_site", DEFAULT_SITE)
    if not cfg.get("downloader"):
        raise BigsdbError("bigsdb_downloader not configured: run '{} setup' "
                          "(or ask the admin for 'setup --shared').".format(PROG))
    secure_token_dir(cfg)
    print("User: {}  |  tokens: {}".format(getpass.getuser(), cfg["token_dir"]))
    backend = downloader_backend(cfg, site, args.timeout)
    fetcher = BigsdbFetcher(site, backend)
    organism = args.organism or ask(
        "Organism your account has access to (e.g. klebsiella, diphtheria)")
    if not organism:
        raise BigsdbError("--organism required.")
    database = resolve_database(fetcher, organism, DbKind.SEQDEF)
    if not args.test_only:
        backend.run_setup(database)

    print("\nTesting access...")
    try:
        backend.get_text("{}/db/{}".format(SITES[site], database))
        cfg["auth"][site] = {"ok": True, "last_check": now_str()}
        print("[OK] {} authentication works.".format(site))
        status = 0
    except BigsdbError as exc:
        cfg["auth"][site] = {"ok": False, "last_check": now_str()}
        print("[FAILED] {}".format(exc))
        print("Also check that your account is registered for the {} database (on the website)."
              .format(database))
        status = 1
    set_user(cfg, "auth", cfg["auth"])
    secure_token_dir(cfg)  # restrict the files bigsdb_downloader has just written
    save_config(cfg_file, cfg)
    return status


def cmd_status(args: argparse.Namespace, cfg: Config, cfg_file: Path) -> int:
    meta = cfg["_meta"]

    def origin(key: str) -> str:
        if key in meta["user"]:
            return "personal"
        return "shared" if key in meta["shared"] else "default"

    print("\n# BIGSip {} status #".format(__version__))
    print("User            : {}".format(getpass.getuser()))
    print("Shared config   : {}{}".format(meta["shared_path"],
                                          "" if meta["shared"] else " (none)"))
    print("Personal config : {}{}".format(cfg_file, "" if cfg_file.exists() else " (none)"))
    print("Downloader      : {}  [{}]".format(cfg.get("downloader") or "not configured",
                                              origin("downloader")))
    if cfg.get("downloader_python"):
        print("Python (.py)    : {}  [{}]".format(cfg["downloader_python"],
                                                  origin("downloader_python")))
    print("Token directory : {}  [{}]{}".format(
        cfg["token_dir"], origin("token_dir"),
        "  (SHARED ACCOUNT)" if cfg.get("shared_tokens") else ""))
    print("Default site    : {}  [{}]".format(cfg.get("default_site"), origin("default_site")))
    for site in SITES:
        auth = cfg["auth"].get(site, {})
        print("  {:<8} key='{}'  auth={}  last check={}".format(
            site, cfg["key_names"].get(site, site), "OK" if auth.get("ok") else "no",
            auth.get("last_check") or "-"))
    return 0


# --------------------------------------------------------------------------
# Commands: listings
# --------------------------------------------------------------------------
def cmd_organisms(args: argparse.Namespace, cfg: Config, cfg_file: Path) -> int:
    fetcher = meta_fetcher(args, cfg)
    need = {"genomes": DbKind.ISOLATES, "loci": DbKind.SEQDEF, "all": None}[args.data]
    print_organisms(fetcher.list_organisms(), args.search, need, fetcher.site)
    return 0


def cmd_schemes(args: argparse.Namespace, cfg: Config, cfg_file: Path) -> int:
    """
    List schemes.

    - no --organism            : scheme IDs of every organism with a seqdef database;
    - --organism X             : schemes of X, with a Profiles column;
    - --organism X --scheme N  : loci of scheme N (--save writes them to a file).
    """
    fetcher = meta_fetcher(args, cfg)
    if args.scheme and not args.organism:
        raise BigsdbError("--scheme requires --organism.")

    if args.organism and args.scheme:
        database = resolve_database(fetcher, args.organism, DbKind.SEQDEF)
        validate_schemes(fetcher, database, organism_of(args.organism), [args.scheme])
        info = fetcher.scheme_info(database, args.scheme)
        print("\nScheme {} - {}: {} loci".format(info.id, info.description, len(info.loci)))
        if args.save:
            Path(args.save).write_text("\n".join(info.loci) + "\n", encoding="utf-8")
            print("List saved to {} (reusable with 'loci --loci {}')."
                  .format(args.save, args.save))
        else:
            for locus in info.loci:
                print("  " + locus)
        return 0

    if args.organism:
        resolve_database(fetcher, args.organism, DbKind.SEQDEF)  # checks the name first
        targets = [organism_of(args.organism)]
    else:
        targets = [o.name for o in filter_organisms(fetcher.list_organisms(), need=DbKind.SEQDEF)]
        print("Schemes of {} organism(s) on {} (--organism X for a single one)...".format(
            len(targets), fetcher.site))
    for org in targets:
        database = resolve_database(fetcher, org, DbKind.SEQDEF)
        try:
            schemes = fetcher.list_schemes_with_profiles(database)
        except Exception as exc:
            print("\n{}: scheme list unavailable ({})".format(org, short_message(str(exc), 80)))
            continue
        if args.search:
            schemes = [s for s in schemes if args.search.lower() in s.description.lower()]
            if not schemes and not args.organism:
                continue
        print_schemes(schemes, org)
    print("\nUse the ID with: loci --organism <org> --scheme <ID>   "
          "or   profiles --organism <org> --scheme <ID>")
    return 0


# --------------------------------------------------------------------------
# Commands: downloads
# --------------------------------------------------------------------------
def require_organism(args: argparse.Namespace, cfg: Config, kind: DbKind) -> None:
    """--organism is mandatory for downloads: if missing, show the list and stop."""
    if args.organism:
        return
    print("\n[ERROR] Missing organism: use --organism <name>.")
    try:
        fetcher = meta_fetcher(args, cfg)
        print_organisms(fetcher.list_organisms(), need=kind, site=fetcher.site)
    except Exception as exc:
        print("(organism list unavailable: {})".format(short_message(str(exc), 100)))
    raise BigsdbError("Missing organism.")


def show_schemes_and_stop(args: argparse.Namespace, cfg: Config, what: str,
                          profiles_only: bool = False) -> None:
    """Called when no scheme was given: show the usable IDs, then stop."""
    fetcher = meta_fetcher(args, cfg)
    database = resolve_database(fetcher, args.organism, DbKind.SEQDEF)
    org = organism_of(args.organism)
    if not (profiles_only and print_profile_schemes(fetcher, database, org) is not None):
        print_schemes(fetcher.list_schemes(database), org)
    raise BigsdbError("Give {} (see the IDs above).".format(what))


def collect_ids(args: argparse.Namespace) -> List[int]:
    tokens: List[str] = []
    if args.ids:
        tokens += names_or_file(args.ids, numeric=True)
    if args.file:
        tokens += read_list_file(args.file)
    if not tokens:
        raise BigsdbError("Give isolates with --ids (e.g. 1,2,10-20 or file.txt) or --file.")
    return expand_ids(tokens)


def cmd_genomes(args: argparse.Namespace, cfg: Config, cfg_file: Path) -> int:
    require_organism(args, cfg, DbKind.ISOLATES)
    ids = collect_ids(args)
    fetcher = make_fetcher(args, cfg)
    ok = do_genomes(fetcher, args.organism, ids, args.outdir, args.header, args.merge,
                    args.overwrite, args.jobs, args.dry_run)
    return 0 if ok else 2


def cmd_loci(args: argparse.Namespace, cfg: Config, cfg_file: Path) -> int:
    """Allele FASTA files: --scheme ID (every locus) and/or --loci a,b,c / file."""
    require_organism(args, cfg, DbKind.SEQDEF)
    loci = names_or_file(args.loci) if args.loci else []
    if args.file:
        loci += read_list_file(args.file)
    schemes = collect_scheme_ids(args.scheme)
    if not loci and not schemes:
        show_schemes_and_stop(args, cfg, "--scheme <ID> (every locus of the scheme) "
                                         "or --loci a,b,c / file")
    fetcher = make_fetcher(args, cfg)
    ok = do_loci(fetcher, args.organism, loci, schemes, args.outdir, args.overwrite,
                 args.jobs, args.dry_run)
    return 0 if ok else 2


def cmd_profiles(args: argparse.Namespace, cfg: Config, cfg_file: Path) -> int:
    require_organism(args, cfg, DbKind.SEQDEF)
    schemes = collect_scheme_ids(args.scheme)
    if not schemes:
        show_schemes_and_stop(args, cfg, "at least one --scheme <ID>", profiles_only=True)
    fetcher = make_fetcher(args, cfg)
    ok = do_profiles(fetcher, args.organism, schemes, args.outdir, args.overwrite, args.dry_run)
    return 0 if ok else 2


# --------------------------------------------------------------------------
# Command: interactive wizard
# --------------------------------------------------------------------------
WIZARD_ACTIONS = (
    ("genomes", "Isolate genomes (contigs FASTA)"),
    ("loci", "Allele FASTA files of loci (whole scheme or list of loci)"),
    ("profiles", "Profiles of a scheme"),
    ("organisms", "Just list the available organisms"),
    ("schemes", "Just list the schemes (IDs) of an organism"),
)


def _wizard_genomes(fetcher: BigsdbFetcher, outdir: str, base_cmd: str, jobs: int) -> int:
    org = pick_organism(fetcher, DbKind.ISOLATES)
    print("\nNUMERIC isolate IDs (not strain names): '1,2,3', ranges '10-20',"
          "\nor path to a text file (one ID per line).")

    def parse(value: str) -> Tuple[List[int], str]:
        if Path(value).expanduser().is_file():
            return expand_ids(read_list_file(value)), "--ids {}".format(value)
        return expand_ids(split_items(value)), "--ids {}".format(",".join(split_items(value)))

    ids, option = ask_until("Isolate ID(s)", parse)
    merge = ask_yes("Also create a merged multi-FASTA?", False)
    print("\nEquivalent command:\n  {} {} --outdir {}{}\n".format(
        base_cmd.format("genomes", org), option, outdir, " --merge" if merge else ""))
    return 0 if do_genomes(fetcher, org, ids, outdir, merge=merge, jobs=jobs) else 2


def _wizard_loci(fetcher: BigsdbFetcher, outdir: str, base_cmd: str, jobs: int) -> int:
    org = pick_organism(fetcher, DbKind.SEQDEF)
    print("\nLoci source: 1. a whole scheme   2. list of loci / file")
    if ask("Choice", "1") == "1":
        database = resolve_database(fetcher, org, DbKind.SEQDEF)
        print_schemes(fetcher.list_schemes(database), org)
        schemes = ask_until("Scheme ID(s) (numeric, e.g. 4)", lambda v: collect_scheme_ids([v]))
        print("\nEquivalent command:\n  {} --scheme {} --outdir {}\n".format(
            base_cmd.format("loci", org), ",".join(schemes), outdir))
        return 0 if do_loci(fetcher, org, scheme_ids=schemes, outdir=outdir, jobs=jobs) else 2

    print("\nLocus names (a locus name is its BIGSdb identifier, e.g. gapA).")
    value = ask("Comma-separated loci, or path to a file (one locus per line)")
    if Path(value).expanduser().is_file():
        loci, option = read_list_file(value), "--loci {}".format(value)
    else:
        loci, option = split_names(value), '--loci "{}"'.format(",".join(split_names(value)))
    print("\nEquivalent command:\n  {} {} --outdir {}\n".format(
        base_cmd.format("loci", org), option, outdir))
    return 0 if do_loci(fetcher, org, loci=loci, outdir=outdir, jobs=jobs) else 2


def _wizard_profiles(fetcher: BigsdbFetcher, outdir: str, base_cmd: str) -> int:
    org = pick_organism(fetcher, DbKind.SEQDEF)
    database = resolve_database(fetcher, org, DbKind.SEQDEF)
    with_pk = print_profile_schemes(fetcher, database, org)   # only schemes with profiles
    if with_pk is None:
        print_schemes(fetcher.list_schemes(database), org)
    elif not with_pk:
        raise BigsdbError("{} has no scheme with profiles. Use 'loci --scheme' for the "
                          "sequences.".format(org))

    def parse(value: str) -> List[str]:
        ids = collect_scheme_ids([value])
        validate_schemes(fetcher, database, org, ids)
        check_profile_schemes(fetcher, database, org, ids)
        return ids

    example = with_pk[0].id if with_pk else "1"
    schemes = ask_until("Scheme ID(s) (numeric, e.g. {})".format(example), parse)
    print("\nEquivalent command:\n  {} --scheme {} --outdir {}\n".format(
        base_cmd.format("profiles", org), ",".join(schemes), outdir))
    return 0 if do_profiles(fetcher, org, schemes, outdir) else 2


def cmd_wizard(args: argparse.Namespace, cfg: Config, cfg_file: Path) -> int:
    """Step-by-step wizard; prints the equivalent command before each download."""
    print("=" * 60)
    print(" BIGSip {} - BIGSdb download wizard".format(__version__))
    print("=" * 60)
    print(POLICY_TEXT + "\n")

    if not cfg.get("downloader") and not args.public:
        print("bigsdb_downloader is not configured yet.")
        if ask_yes("Configure it now?"):
            setup_args = argparse.Namespace(downloader=None, token_dir=None, site=None,
                                            key_name=None, shared=False, shared_tokens=False,
                                            python=None)
            cmd_setup(setup_args, cfg, cfg_file)
            cfg = load_config(cfg_file)
        else:
            print("-> Public mode (pre-2025 data only).")
            args.public = True

    site = args.site or ask("Site (Pasteur/PubMLST)", cfg.get("default_site", DEFAULT_SITE))
    if site not in SITES:
        raise BigsdbError("Unknown site: {}".format(site))
    args.site = site

    if not args.public and not cfg["auth"].get(site, {}).get("ok"):
        print("\n{} authentication has not been validated.".format(site))
        if ask_yes("Run the OAuth setup now?"):
            cmd_auth(argparse.Namespace(site=site, organism=None, timeout=args.timeout,
                                        test_only=False), cfg, cfg_file)
            cfg = load_config(cfg_file)

    print("\nWhat do you want to download?")
    action = choose_from(WIZARD_ACTIONS, lambda a: a[1], "Number", key=lambda a: a[0])[0]
    fetcher = make_fetcher(args, cfg)

    if action == "organisms":
        print_organisms(fetcher.list_organisms(), ask("Filter (empty = show all)"), site=site)
        return 0
    if action == "schemes":
        org = pick_organism(fetcher, DbKind.SEQDEF)
        database = resolve_database(fetcher, org, DbKind.SEQDEF)
        print_schemes(fetcher.list_schemes_with_profiles(database), org)
        return 0

    outdir = ask("Output directory", args.outdir)
    base_cmd = "python {} {{}} --site {} --organism {{}}{}".format(
        Path(sys.argv[0]).name, site, " --public" if args.public else "")
    if action == "genomes":
        return _wizard_genomes(fetcher, outdir, base_cmd, args.jobs)
    if action == "loci":
        return _wizard_loci(fetcher, outdir, base_cmd, args.jobs)
    return _wizard_profiles(fetcher, outdir, base_cmd)


# --------------------------------------------------------------------------
# Command-line interface
# --------------------------------------------------------------------------
EX_START = """
  Configuration
    {p} setup
    {p} auth --organism diphtheria
""".format(p=PROG)

EX_FIND = """
  Organisms and schemes
    {p} organisms                                # organisms
    {p} schemes                                  # schemes, all
    {p} schemes --organism klebsiella            # schemes of klebsiella
    {p} schemes --organism klebsiella --scheme 4 # loci of scheme 4
""".format(p=PROG)

EX_GENOMES = """
  Genomes
    {p} genomes --organism klebsiella --ids 1,2,10-20
    {p} genomes --organism klebsiella --ids ids.txt
""".format(p=PROG)

EX_LOCI = """
  Loci (alleles)
    {p} loci --organism klebsiella --scheme 4
    {p} loci --organism klebsiella --loci gapA,infB
    {p} loci --organism klebsiella --loci loci.txt
""".format(p=PROG)

EX_PROFILES = """
  Profiles
    {p} profiles --organism klebsiella --scheme 4
""".format(p=PROG)

EX_SHARED = """
  Shared install (admin)
    {p} setup --shared --downloader /opt/venv/bin/bigsdb-downloader
""".format(p=PROG)

RAW = argparse.RawDescriptionHelpFormatter


def build_parser() -> argparse.ArgumentParser:
    common = argparse.ArgumentParser(add_help=False)
    group = common.add_argument_group("common options")
    group.add_argument("--site", choices=list(SITES),
                       help="BIGSdb site (default: the one in the config)")
    group.add_argument("--config",
                       help="Config file (default: {} or ${})".format(DEFAULT_CONFIG, CONFIG_ENV))
    group.add_argument("--public", action="store_true",
                       help="No authentication (data produced before 01/01/2025 only)")
    group.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT,
                       help="Timeout per request in s (default {})".format(DEFAULT_TIMEOUT))
    verbosity = group.add_mutually_exclusive_group()
    verbosity.add_argument("-v", "--verbose", action="store_true",
                           help="Show debug messages (commands run, fallbacks)")
    verbosity.add_argument("-q", "--quiet", action="store_true",
                           help="Only show warnings, errors and summaries")

    download = argparse.ArgumentParser(add_help=False)
    group = download.add_argument_group("download options")
    group.add_argument("--organism", "--db", dest="organism", metavar="ORGANISM",
                       help="Organism name (e.g. klebsiella, bordetella) - see 'organisms'")
    group.add_argument("--outdir", default=DEFAULT_OUTDIR,
                       help="Output directory (default: ./%(default)s)")
    group.add_argument("--overwrite", action="store_true",
                       help="Download again files that are already present")
    group.add_argument("--dry-run", action="store_true",
                       help="Show the commands without running them")

    parser = argparse.ArgumentParser(
        prog="bigsip.py",
        description="BIGSip - sip just the BIGSdb data you need, through bigsdb_downloader.\n\n"
                    + POLICY_TEXT,
        epilog="Examples:\n" + EX_START + EX_FIND + EX_GENOMES + EX_LOCI + EX_PROFILES,
        formatter_class=RAW, parents=[common])
    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)
    sub = parser.add_subparsers(dest="command", metavar="command")

    p = sub.add_parser("setup", parents=[common], formatter_class=RAW,
                       help="Set the bigsdb_downloader path and the token directory",
                       epilog="Examples:" + EX_START + EX_SHARED)
    p.add_argument("--downloader",
                   help="Path to bigsdb_downloader.py or to the bigsdb-downloader executable")
    p.add_argument("--token-dir", help="Directory where keys/tokens are stored")
    p.add_argument("--key-name", help="API key name for the site (default: the site name)")
    p.add_argument("--python", help="Python to use when --downloader is a .py file "
                                    "(must have 'rauth'; e.g. /opt/venv/bin/python)")
    p.add_argument("--shared", action="store_true",
                   help="ADMIN: write the shared config next to the script (downloader + token "
                        "directory template, e.g. ~/.bigsip/tokens or /data/tokens/{user})")
    p.add_argument("--shared-tokens", action="store_true",
                   help="ADMIN: allow a common token directory (single service account)")
    p.set_defaults(func=cmd_setup)

    p = sub.add_parser("auth", parents=[common],
                       help="Run the OAuth authentication (interactive) and test it")
    p.add_argument("--organism", "--db", dest="organism", metavar="ORGANISM",
                   help="Organism your account has access to (needed for the setup)")
    p.add_argument("--test-only", action="store_true",
                   help="Only test, do not run the setup again")
    p.set_defaults(func=cmd_auth)

    p = sub.add_parser("status", parents=[common], help="Show the configuration")
    p.set_defaults(func=cmd_status)

    p = sub.add_parser("organisms", aliases=["databases"], parents=[common], formatter_class=RAW,
                       help="List the organisms available on a site",
                       epilog="Examples:" + EX_FIND)
    p.add_argument("--search", help="Text filter (name or description)")
    p.add_argument("--data", choices=["genomes", "loci", "all"], default="all",
                   help="Only show organisms that have this type of data")
    p.set_defaults(func=cmd_organisms)

    p = sub.add_parser("schemes", parents=[common], formatter_class=RAW,
                       help="List scheme IDs (all organisms, or a single one)",
                       epilog="Examples:" + EX_FIND)
    p.add_argument("--organism", "--db", dest="organism", metavar="ORGANISM",
                   help="One organism (default: all)")
    p.add_argument("--search", help="Filter on the scheme description (e.g. cgMLST, MLST)")
    p.add_argument("--scheme", help="Show the loci of this scheme")
    p.add_argument("--save", help="With --scheme: save the list of loci to a file")
    p.set_defaults(func=cmd_schemes)

    p = sub.add_parser("genomes", parents=[common, download], formatter_class=RAW,
                       help="Download isolate genomes", epilog="Examples:" + EX_GENOMES)
    p.add_argument("--ids", help="NUMERIC isolate IDs (not strain names): 1,2,10-20, "
                                 "or a file (one ID per line)")
    p.add_argument("--file", help="Text file, one ID per line")
    p.add_argument("--header", choices=["id", "original_designation"],
                   help="Contig headers in the FASTA")
    p.add_argument("--merge", action="store_true", help="Also create a single multi-FASTA")
    p.add_argument("--jobs", type=int, default=1,
                   help="Parallel downloads (default 1, 4 max recommended)")
    p.set_defaults(func=cmd_genomes)

    p = sub.add_parser("loci", aliases=["alleles"], parents=[common, download],
                       formatter_class=RAW, epilog="Examples:" + EX_LOCI,
                       help="Download allele FASTA files: of a whole scheme or of chosen loci")
    p.add_argument("--scheme", action="append",
                   help="NUMERIC scheme ID(s): downloads EVERY locus of the scheme "
                        "(3,4 or repeatable)")
    p.add_argument("--loci", help="Locus names (name = BIGSdb identifier), comma-separated, "
                                  "OR a text file (one locus per line)")
    p.add_argument("--file", help="Text file, one locus per line (same as --loci file.txt)")
    p.add_argument("--jobs", type=int, default=1, help="Parallel downloads (default 1)")
    p.set_defaults(func=cmd_loci)

    p = sub.add_parser("profiles", parents=[common, download], formatter_class=RAW,
                       help="Download scheme profiles", epilog="Examples:" + EX_PROFILES)
    p.add_argument("--scheme", action="append",
                   help="NUMERIC scheme ID(s) (3,4 or repeatable); "
                        "without --scheme: shows the available IDs")
    p.set_defaults(func=cmd_profiles)

    p = sub.add_parser("wizard", parents=[common],
                       help="Interactive wizard (default when no argument)")
    p.add_argument("--outdir", default=DEFAULT_OUTDIR)
    p.add_argument("--jobs", type=int, default=1)
    p.set_defaults(func=cmd_wizard)
    return parser


def setup_logging(verbose: bool = False, quiet: bool = False) -> None:
    """Messages go to stderr as '[LEVEL] message'; tables and summaries stay on stdout."""
    level = logging.DEBUG if verbose else logging.WARNING if quiet else logging.INFO
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
    log.handlers[:] = [handler]
    log.setLevel(level)
    log.propagate = False


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    argv = list(sys.argv[1:] if argv is None else argv)
    args = parser.parse_args(argv)
    if not args.command:
        args = parser.parse_args(["wizard"] + argv)
    for name, default in (("dry_run", False), ("jobs", 1), ("outdir", DEFAULT_OUTDIR),
                          ("python", None)):
        if not hasattr(args, name):
            setattr(args, name, default)
    if args.jobs < 1:
        parser.error("--jobs must be >= 1")
    setup_logging(args.verbose, args.quiet)

    cfg_file = config_path(args.config)
    cfg = load_config(cfg_file)
    try:
        return args.func(args, cfg, cfg_file)
    except BigsdbError as exc:
        log.error("%s", exc)
        return 1
    except KeyboardInterrupt:
        log.error("Interrupted. Files already downloaded are kept; "
                  "run the same command again to resume.")
        return 130


if __name__ == "__main__":
    sys.exit(main())