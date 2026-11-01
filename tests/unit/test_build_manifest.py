# tests/unit/test_build_manifest.py
import hashlib
import json
import subprocess
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

from chd_atlas.build.emit import Emitter
from chd_atlas.build.manifest import SCHEMA_VERSION, STATUS, source_commit, write_manifest
from chd_atlas.corpus import Corpus


def _published(root: Path) -> dict[str, Any]:
    """The manifest as a consumer gets it: parsed from the file that was written."""
    manifest: dict[str, Any] = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    return manifest


def _git(root: Path, *args: str) -> str:
    completed = subprocess.run(
        ["git", "-C", str(root), *args], capture_output=True, text=True, check=True
    )
    return completed.stdout.strip()


def _checkout(root: Path, *, commits: int) -> Path:
    """A real repository, built without depending on the machine's git identity.

    The `-c` overrides matter: a machine configured with `commit.gpgsign = true`
    and no key, or with no `user.email` at all, would otherwise fail to commit
    and turn a fixture into a red test about something else entirely.
    """
    root.mkdir(parents=True, exist_ok=True)
    _git(root, "init", "-q")
    for number in range(commits):
        _git(
            root,
            "-c",
            "user.name=curator",
            "-c",
            "user.email=curator@example.com",
            "-c",
            "commit.gpgsign=false",
            "commit",
            "-q",
            "--allow-empty",
            "-m",
            f"commit {number}",
        )
    return root


def test_the_manifest_advertises_the_digest_of_the_bytes_on_disk(tmp_path: Path) -> None:
    """The property the manifest exists for, and the one easiest to fake passing.

    A checksum that disagrees with the file is worse than no manifest: a
    consumer verifying its download rejects a file that is correct. The
    specified suite asserted only that each value began with `sha256:` and
    passed with every checksum replaced by 64 zeros (A43), so this recomputes
    the digest from the bytes that were actually written.

    `hashlib` directly rather than `emit.checksum`: a suite that hashes with the
    same function it is checking cannot see the two disagree.

    A gzipped artifact alongside a plain one because that is the pair that can
    diverge — for a `.gz` the served bytes are the compressed ones, and
    checksumming the JSON inside would publish a value no consumer can verify.

    The set equality is also what pins the self-exclusion: `manifest.json`
    cannot appear in its own `files`, since the value would have to be computed
    from bytes containing it.
    """
    emitter = Emitter(root=tmp_path)
    emitter.write_json("genes/index.json", {"genes": []})
    emitter.write_json_gz("variants/1.json.gz", {"variants": []})

    write_manifest(Corpus(root=Path(".")), emitter, commit="abc123", counts={})

    files = _published(tmp_path)["files"]
    assert set(files) == {"genes/index.json", "variants/1.json.gz"}
    for relative, advertised in files.items():
        served = (tmp_path / relative).read_bytes()
        assert advertised == "sha256:" + hashlib.sha256(served).hexdigest(), relative


def test_the_manifest_counts_every_record_kind_separately(tmp_path: Path) -> None:
    """Nine kinds, nine different lengths, so a cross-wired count cannot pass.

    `write_manifest` only ever calls `len()` over the corpus, so the fixture
    varies lengths rather than record content; the placeholders stand in for
    records because nothing here reads a field of one. With equal lengths
    `"assertions": len(corpus.datasets)` would publish the right number and the
    swap would be invisible.

    The three build counts (2.7) are asserted in the same equality rather than
    separately, because the property that matters is the *whole* census: a
    consumer summarising the atlas reads this object, and through 2.6 it
    published `assertions: 1, datasets: 0, functional: 0` for a site publishing
    23 genes and 290 burden statistics. They are given values distinct from every
    corpus length here for the same reason the corpus lengths differ from each
    other.

    The lengths start at two rather than one for the reason A30 records: `1 ==
    True` in Python, so a count of one is also satisfied by a `bool()` over the
    same tuple, and whole-dict equality would accept it.
    """
    corpus = Corpus(  # type: ignore[arg-type]
        root=Path("."),
        assertions=("a",) * 2,
        datasets=("d",) * 3,
        featured=("f",) * 4,
        functional=("u",) * 5,
        phenotypes=("p",) * 6,
        publications=("b",) * 7,
        # 11, not 10: `cohorts` (the registry `cohorts.json` publishes) and
        # `cohort_families` (how many independent collections the burden rows
        # group into) are different quantities with confusable names, and giving
        # them one value would let each be published under the other's name
        # undetected -- CLAUDE.md section 4.30.
        cohorts=("c",) * 11,
        # 12: the scope registry `scope.json` publishes (schema 2.15), distinct
        # from every other length so a cross-wire cannot pass.
        chd_scope=("s",) * 12,
    )
    emitter = Emitter(root=tmp_path)

    write_manifest(
        corpus,
        emitter,
        commit="abc123",
        counts={"genes": 8, "burden_rows": 9, "cohort_families": 10},
    )

    assert _published(tmp_path)["counts"] == {
        "assertions": 2,
        "cohorts": 11,
        "datasets": 3,
        "featured": 4,
        "functional": 5,
        "phenotypes": 6,
        "publications": 7,
        "scope_terms": 12,
        "genes": 8,
        "burden_rows": 9,
        "cohort_families": 10,
    }


def test_a_build_count_may_not_shadow_a_corpus_count(tmp_path: Path) -> None:
    """`counts` is merged, so a colliding key would silently republish one figure.

    `datasets` is the live hazard rather than a hypothetical one: it means
    *omics* datasets and is 0, and the burden census's third figure is a count of
    independent datasets. A caller reaching for the obvious name would overwrite
    a corpus figure with a build one, and `sort_keys` would leave the output
    looking entirely ordinary — no test reading a single key could see it.

    Refused rather than namespaced, because renaming the incoming key is the fix
    and a nested object would make every consumer walk two levels for one census.
    """
    emitter = Emitter(root=tmp_path)

    with pytest.raises(ValueError, match="would shadow the corpus counts"):
        write_manifest(Corpus(root=Path(".")), emitter, commit="abc123", counts={"datasets": 3})


@pytest.mark.parametrize("commit", ["a" * 40, None], ids=["a checkout", "a tarball"])
def test_the_manifest_publishes_five_keys_and_nothing_that_varies(
    tmp_path: Path, commit: str | None
) -> None:
    """No build timestamp, by decision, and an empty corpus still gets a manifest.

    A timestamp would make two builds of one commit differ, so no test could
    assert byte-identical rebuilds and the checksums would stop being a pure
    function of the inputs; a consumer wanting a date reads the commit's. Pinning
    the exact key set is what makes adding one a decision rather than a side
    effect.

    `None` is a row rather than a test of its own because it is the tarball build
    that `source_commit` returns None for, and what matters is that it reaches
    the file as JSON `null` — a key a consumer can test — rather than being
    omitted or stringified to "None" on the way in.

    `schema_version` is asserted against a literal as well as the constant: the
    constant compared with itself passes whatever it is set to, and the version
    is a published contract a consumer branches on. `status` is asserted the
    same way, against `STATUS` and against the literal `"in-development"` a
    consumer actually branches on — a programmatic one now, not only a reader of
    `index.html`'s prose, which is the field's whole reason to exist.
    """
    emitter = Emitter(root=tmp_path)

    write_manifest(Corpus(root=Path(".")), emitter, commit=commit, counts={})

    manifest = _published(tmp_path)
    assert set(manifest) == {"schema_version", "source_commit", "status", "counts", "files"}
    assert manifest["schema_version"] == "2.20" == SCHEMA_VERSION
    assert manifest["status"] == "in-development" == STATUS
    assert manifest["source_commit"] == commit
    assert manifest["files"] == {}


def test_nothing_can_be_written_once_the_manifest_is_out(tmp_path: Path) -> None:
    """The manifest is written last — enforced here rather than merely intended.

    A file written afterwards is on disk and served while appearing in no
    manifest and carrying no checksum: invisible to every consumer that verifies
    what it downloaded, and afterwards indistinguishable from a file that
    belonged there. Reproduced before the guard existed — `late/shard.json` was
    written, served, and nothing raised (A44).

    A second manifest fails through the same guard, which is the right answer for
    the same reason: it would advertise a checksum for the first, and only one of
    the two is the file a consumer fetched.
    """
    emitter = Emitter(root=tmp_path)
    emitter.write_json("genes/index.json", {"genes": []})
    write_manifest(Corpus(root=Path(".")), emitter, commit="abc123", counts={})
    published = dict(emitter.checksums)

    with pytest.raises(ValueError, match="after the manifest"):
        emitter.write_json("late/shard.json", {"late": True})
    with pytest.raises(ValueError, match="after the manifest"):
        write_manifest(Corpus(root=Path(".")), emitter, commit="abc123", counts={})

    # Refused has to mean not written, not written-and-forgotten: a file on disk
    # that the manifest does not list is the whole failure being prevented.
    assert not (tmp_path / "late").exists()
    assert emitter.checksums == published


def test_source_commit_reads_the_head_of_the_checkout_it_is_given(tmp_path: Path) -> None:
    """Provenance is the commit at `root`, not wherever the build process runs.

    The atlas's own suite runs inside a checkout, so a `source_commit` that
    ignored `root` would answer with *this* repository's HEAD and look entirely
    plausible. A second repository one directory over is what makes that visible,
    and is the reason `-C` is passed at all.

    Compared against the full sha read back from git rather than against a shape,
    which is what separates `rev-parse HEAD` from `rev-parse --short HEAD` and
    pins that the newline git prints is stripped. The length assertion keeps that
    honest without assuming sha1: a short sha is 7 or 8 characters, sha1 is 40
    and sha256 is 64, so the comparison is only meaningful if it is at least 40.
    """
    checkout = _checkout(tmp_path / "elsewhere", commits=2)
    head = _git(checkout, "rev-parse", "HEAD")

    assert len(head) >= 40
    assert source_commit(checkout) == head


def test_a_checkout_whose_path_contains_a_space_still_reports_its_commit(
    tmp_path: Path,
) -> None:
    """`--show-toplevel` prints a path, and paths have spaces in them.

    `~/My Projects/atlas` and `~/My Drive/…` are ordinary on macOS. Parsing the
    two-line output by whitespace made such a path three tokens instead of two,
    so a perfectly valid checkout published no provenance — and did it by
    failing into the same `None` that means "not a checkout at all", which is
    the reading a curator would take from `"source_commit": null`.

    Asserted against the real HEAD rather than against not-None, because the
    line-based parse has to put the two fields the right way round: returning the
    *path* as the commit would also be non-null.
    """
    checkout = _checkout(tmp_path / "My Projects" / "atlas", commits=1)

    assert " " in str(checkout)
    assert source_commit(checkout) == _git(checkout, "rev-parse", "HEAD")


def test_a_root_inside_an_unrelated_checkout_has_no_provenance(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """`git -C` searches upward, so "not a checkout" has to be made to mean it.

    A tarball unpacked anywhere inside some other repository — a scratch clone, a
    monorepo, a home directory under version control — is not a checkout of this
    atlas, and the docstring and `docs/data-api.md` both promise `null` for it.
    Before this guard it published the *enclosing* repository's HEAD: a real,
    verifiable sha naming a commit that contains none of the built data.

    That is strictly worse than the `"HEAD"` the returncode check exists to
    prevent, and for the opposite reason. A provenance claim resolving to nothing
    announces itself the first time anyone follows it; one that resolves to a
    real commit in the wrong repository never does.

    `GIT_DIR` is set here because it overrides `-C` outright, so a build
    inheriting it from a surrounding tool would report that repository's commit
    however carefully `root` was chosen. It is the same defect by a quieter
    route, which is why the environment is cleared rather than trusted.
    """
    enclosing = _checkout(tmp_path / "enclosing", commits=1)
    unpacked = enclosing / "tarball"
    unpacked.mkdir()

    assert source_commit(unpacked) is None
    # The enclosing repository is genuinely resolvable — the None above is the
    # guard working, not git failing to find anything.
    assert source_commit(enclosing) == _git(enclosing, "rev-parse", "HEAD")

    monkeypatch.setenv("GIT_DIR", str(enclosing / ".git"))
    assert source_commit(unpacked) is None


def _plain_directory(tmp_path: Path) -> Path:
    return tmp_path


def _repository_with_no_commits(tmp_path: Path) -> Path:
    return _checkout(tmp_path / "fresh", commits=0)


def _directory_that_does_not_exist(tmp_path: Path) -> Path:
    return tmp_path / "gone"


@pytest.mark.parametrize(
    "build_root",
    [_plain_directory, _repository_with_no_commits, _directory_that_does_not_exist],
    ids=["not a checkout", "checkout with no commits", "root does not exist"],
)
def test_a_root_with_no_resolvable_head_has_no_provenance(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, build_root: Callable[[Path], Path]
) -> None:
    """A tarball build must still produce a complete site, just an unprovenanced one.

    The second row is not a variation on the first. In a repository with no
    commits `git rev-parse HEAD` exits 128 *and prints the literal string `HEAD`
    on stdout*, so a `source_commit` that trusted stdout alone would publish
    `"source_commit": "HEAD"` — provenance that resolves to nothing, which a
    consumer cannot tell apart from a real sha. The third row pins that pointing
    at a directory that is not there returns rather than raises.
    """
    # `git rev-parse` walks *up* the tree, so "this directory is not in a
    # checkout" is a fact about the machine unless the test forces it: with
    # $TMPDIR inside a repository the first row would report that repository's
    # HEAD. Measured here — pytest's tmp_path is under /private/var and is in no
    # checkout, so the row passes either way today — but the ceiling makes the
    # precondition part of the test instead of part of the environment. It stops
    # the walk above tmp_path, so the second row's repository, which is below it,
    # is still found.
    monkeypatch.setenv("GIT_CEILING_DIRECTORIES", str(tmp_path.parent))

    assert source_commit(build_root(tmp_path)) is None


@pytest.mark.parametrize(
    "answer",
    [
        FileNotFoundError(2, "No such file or directory: 'git'"),
        subprocess.TimeoutExpired(["git"], 10),
        subprocess.CompletedProcess(["git"], returncode=0, stdout="\n", stderr=""),
    ],
    ids=["git is not installed", "git hung", "git answered with nothing"],
)
def test_a_git_that_does_not_answer_leaves_the_build_unprovenanced(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, answer: object
) -> None:
    """Three failures no repository fixture can stage, and one of them is silent.

    A missing git raises `FileNotFoundError` from `subprocess.run` rather than
    returning non-zero, and a machine without git must still be able to publish —
    the same decision as the tarball build.

    An empty answer with a zero exit is the one shape the returncode check cannot
    cover, and it must become `None` rather than `""`: an empty string is falsy
    but still reads as a recorded commit to a consumer testing the key for null.

    The timeout is asserted because nothing else can pin it — the hung row proves
    only that the exception is caught, not that anything raises it — and dropping
    it turns a git blocked on a dead network mount or a credential prompt into a
    build that never finishes and never says why.
    """
    calls: list[dict[str, Any]] = []

    def _fake_run(command: list[str], **kwargs: Any) -> subprocess.CompletedProcess[str]:
        calls.append(kwargs)
        if isinstance(answer, BaseException):
            raise answer
        assert isinstance(answer, subprocess.CompletedProcess)
        return answer

    monkeypatch.setattr(subprocess, "run", _fake_run)

    assert source_commit(tmp_path) is None
    assert calls[0]["timeout"] == 10


def test_the_manifest_states_the_cell_type_census(tmp_path: Path) -> None:
    """`cell_type_datasets` and `cell_type_genes` are published as given, distinct."""
    emitter = Emitter(root=tmp_path)

    write_manifest(
        Corpus(root=Path(".")),
        emitter,
        commit=None,
        counts={"cell_type_genes": 7, "cell_type_datasets": 3},
    )

    counts = _published(tmp_path)["counts"]
    assert counts["cell_type_genes"] == 7
    assert counts["cell_type_datasets"] == 3
