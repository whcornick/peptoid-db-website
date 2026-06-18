import re


CODE_PATTERN = re.compile(
    r"^(?P<year>\d{2})(?P<publication>[A-Z]{2})"
    r"(?P<structure>\d+)-(?P<length>\d+)-(?P<topology>[ACM])$"
)


class CodeGenerationError(Exception):
    pass


def _normalize_doi(value):
    value = (value or "").strip().lower()
    for prefix in (
        "https://doi.org/",
        "http://doi.org/",
        "doi:",
    ):
        if value.startswith(prefix):
            value = value[len(prefix):]
    return value.strip()


def _letters_to_index(letters):
    return (
        (ord(letters[0]) - ord("A")) * 26
        + ord(letters[1]) - ord("A")
    )


def _index_to_letters(index):
    if index < 0 or index >= 26 * 26:
        raise CodeGenerationError(
            "No publication identifiers remain for this year."
        )
    return "{}{}".format(
        chr(ord("A") + index // 26),
        chr(ord("A") + index % 26),
    )


def generate_peptoid_code(
    submission, residue_count, peptoids, reserved_submissions=()
):
    doi = _normalize_doi(submission.pub_doi)
    if not doi:
        raise CodeGenerationError(
            "A publication DOI is required for automatic code generation."
        )
    if not submission.release:
        raise CodeGenerationError(
            "A release date is required for automatic code generation."
        )
    if residue_count < 2:
        raise CodeGenerationError(
            "A finalized peptoid must contain at least two residues."
        )

    parsed = []
    matching_publication = []

    for peptoid in peptoids:
        match = CODE_PATTERN.fullmatch(peptoid.code or "")
        if not match:
            continue
        data = match.groupdict()
        parsed.append(data)

        if _normalize_doi(peptoid.pub_doi) == doi:
            matching_publication.append(data)

    reserved = []
    for reserved_submission in reserved_submissions:
        match = CODE_PATTERN.fullmatch(
            reserved_submission.proposed_code or ""
        )
        if not match:
            continue
        data = match.groupdict()
        reserved.append(data)
        parsed.append(data)

        if _normalize_doi(reserved_submission.pub_doi) == doi:
            matching_publication.append(data)

    if matching_publication:
        prefixes = {
            (item["year"], item["publication"])
            for item in matching_publication
        }
        if len(prefixes) != 1:
            raise CodeGenerationError(
                "Existing records for this DOI use inconsistent code prefixes."
            )
        year, publication = prefixes.pop()
    else:
        year = str(submission.release.year)[-2:]
        used = {
            item["publication"]
            for item in parsed
            if item["year"] == year
        }
        next_index = (
            max((_letters_to_index(value) for value in used), default=-1)
            + 1
        )
        publication = _index_to_letters(next_index)

    prefix = "{}{}".format(year, publication)
    used_structure_numbers = [
        int(item["structure"])
        for item in parsed
        if "{}{}".format(item["year"], item["publication"]) == prefix
    ]
    structure_number = max(used_structure_numbers, default=0) + 1

    topology = {
        "linear": "A",
        "cyclic": "C",
        "multicyclic": "M",
        "A": "A",
        "C": "C",
        "M": "M",
    }.get(submission.topology)

    if topology is None:
        raise CodeGenerationError(
            "Unsupported topology: {}".format(submission.topology)
        )

    return "{}{}-{}-{}".format(
        prefix,
        structure_number,
        residue_count,
        topology,
    )
