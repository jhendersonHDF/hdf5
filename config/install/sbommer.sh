#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: sbommer.sh [INPUT_ARCHIVE] [OUTPUT_SBOM]

Generate an SPDX 2.3 JSON SBOM for the HDF5 source archive.

Defaults:
  INPUT_ARCHIVE   hdf5-2.1.1.tar.gz
  OUTPUT_SBOM     hdf5-2.1.1.tar.gz.spdx.json

Optional environment variables:
  SBOM_CREATOR_PERSON      Override the "Person:" creator entry
  SBOM_CREATOR_TOOL        Override the "Tool:" creator entry
  SBOM_DATA_LICENSE        Override the SPDX data license
  SBOM_DOCUMENT_NAMESPACE  Override the document namespace URI
  SBOM_DOWNLOAD_LOCATION   Override the package downloadLocation
  SBOM_LICENSE_LIST_VERSION  Override the SPDX license list version
  SBOM_SPDX_VERSION        Override the SPDX version
EOF
}

require_cmd() {
  local cmd=$1
  if ! command -v "$cmd" >/dev/null 2>&1; then
    printf 'Missing required command: %s\n' "$cmd" >&2
    exit 1
  fi
}

spdx_safe_id_suffix() {
  printf '%s' "$1" | sed -E 's/[^A-Za-z0-9.-]+/-/g; s/^-+//; s/-+$//'
}

archive_stem() {
  local archive_name
  archive_name=$(basename "$1")
  case "$archive_name" in
    *.tar.gz) printf '%s\n' "${archive_name%.tar.gz}" ;;
    *.tgz) printf '%s\n' "${archive_name%.tgz}" ;;
    *.zip) printf '%s\n' "${archive_name%.zip}" ;;
    *) printf '%s\n' "$archive_name" ;;
  esac
}

list_archive_entries() {
  case "$archive_format" in
    tar.gz|tgz) tar -tzf "$input_archive" ;;
    zip) unzip -Z1 "$input_archive" ;;
  esac
}

extract_archive() {
  case "$archive_format" in
    tar.gz|tgz) tar -xzf "$input_archive" -C "$tmpdir" ;;
    zip) unzip -q "$input_archive" -d "$tmpdir" ;;
  esac
}

print_archive_file() {
  local member_path=$1
  case "$archive_format" in
    tar.gz|tgz)
      tar -xOf "$input_archive" "./$member_path" 2>/dev/null ||
        tar -xOf "$input_archive" "$member_path" 2>/dev/null
      ;;
    zip)
      unzip -p "$input_archive" "$member_path"
      ;;
  esac
}

if [[ ${1:-} == "-h" || ${1:-} == "--help" ]]; then
  usage
  exit 0
fi

for cmd in awk date find id jq mktemp sed sha1sum sha256sum sort; do
  require_cmd "$cmd"
done

input_archive=${1:-hdf5-2.1.1.tar.gz}
output_sbom=${2:-${input_archive##*/}.spdx.json}

if [[ ! -f "$input_archive" ]]; then
  printf 'Input archive not found: %s\n' "$input_archive" >&2
  exit 1
fi

case "$input_archive" in
  *.tar.gz) archive_format=tar.gz ;;
  *.tgz) archive_format=tgz ;;
  *.zip) archive_format=zip ;;
  *)
    printf 'Unsupported archive format: %s\n' "$input_archive" >&2
    exit 1
    ;;
esac

case "$archive_format" in
  tar.gz|tgz) require_cmd tar ;;
  zip) require_cmd unzip ;;
esac

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

root_dir=$(
  list_archive_entries | awk '
    {
      entry = $0
      gsub(/^\.\//, "", entry)
      gsub(/\/$/, "", entry)
      if (entry == "") {
        next
      }
      if (root == "") {
        split(entry, parts, "/")
        root = parts[1]
      }
    }
    END {
      print root
    }'
)
if [[ -z "$root_dir" ]]; then
  printf 'Unable to determine archive root directory from %s\n' "$input_archive" >&2
  exit 1
fi

extract_archive
source_dir=$tmpdir/$root_dir
if [[ ! -d "$source_dir" ]]; then
  printf 'Extracted source directory not found: %s\n' "$source_dir" >&2
  exit 1
fi

spdx_version=${SBOM_SPDX_VERSION:-SPDX-2.3}
data_license=${SBOM_DATA_LICENSE:-CC0-1.0}
license_list_version=${SBOM_LICENSE_LIST_VERSION:-3.22}

archive_name=$(basename "$input_archive")
output_name=$(basename "$output_sbom")
package_name="HDF5"
package_spdx_id="SPDXRef-PACKAGE-hdf5"
document_name="HDF5 SBOM"
version_info=$(
  print_archive_file "$root_dir/README.md" 2>/dev/null |
    awk '/^HDF5 version / { version = $3 } END { print version }'
)
if [[ -z "$version_info" ]]; then
  version_info=$(archive_stem "$archive_name")
  version_info=${version_info#hdf5-}
fi

creator_person=${SBOM_CREATOR_PERSON:-$(id -un)}
creator_tool=${SBOM_CREATOR_TOOL:-sbommer.sh}
created_at=$(date -u +%Y-%m-%dT%H:%M:%SZ)
document_namespace=${SBOM_DOCUMENT_NAMESPACE:-https://www.hdfgroup.org/HDF5/$output_name}
download_location=${SBOM_DOWNLOAD_LOCATION:-https://support.hdfgroup.org/downloads/index.html}
package_sha256=$(sha256sum "$input_archive" | awk '{ print $1 }')
package_purl="pkg:generic/org.hdfgroup/hdf5@$version_info"

files_jsonl=$tmpdir/files.jsonl
relationships_jsonl=$tmpdir/relationships.jsonl
sha1_list=$tmpdir/file-sha1s.txt
package_json=$tmpdir/package.json

touch "$files_jsonl" "$relationships_jsonl" "$sha1_list"

# SPDX package verification code is the SHA1 of the sorted file SHA1 values.
while IFS= read -r -d '' relpath; do
  file_sha1=$(sha1sum "$source_dir/$relpath" | awk '{ print $1 }')
  file_sha256=$(sha256sum "$source_dir/$relpath" | awk '{ print $1 }')
  file_spdx_id="SPDXRef-FILE-$(spdx_safe_id_suffix "$relpath")"

  jq -cn \
    --arg spdxid "$file_spdx_id" \
    --arg sha1 "$file_sha1" \
    --arg sha256 "$file_sha256" \
    --arg filename "$relpath" \
    '{
      SPDXID: $spdxid,
      checksums: [
        {algorithm: "SHA1", checksumValue: $sha1},
        {algorithm: "SHA256", checksumValue: $sha256}
      ],
      fileName: $filename
    }' >> "$files_jsonl"

  jq -cn \
    --arg spdxElementId "$package_spdx_id" \
    --arg relatedSpdxElement "$file_spdx_id" \
    '{
      relatedSpdxElement: $relatedSpdxElement,
      relationshipType: "CONTAINS",
      spdxElementId: $spdxElementId
    }' >> "$relationships_jsonl"

  printf '%s\n' "$file_sha1" >> "$sha1_list"
done < <(cd "$source_dir" && find . -type f -printf '%P\0' | LC_ALL=C sort -z)

package_verification_code=$(LC_ALL=C sort "$sha1_list" | tr -d '\n' | sha1sum | awk '{ print $1 }')

jq -cn \
  --arg spdxid "$package_spdx_id" \
  --arg checksum "$package_sha256" \
  --arg download_location "$download_location" \
  --arg purl "$package_purl" \
  --arg license "BSD-3-Clause" \
  --arg name "$package_name" \
  --arg originator "Organization: The HDF Group" \
  --arg package_file_name "$archive_name" \
  --arg verification_code "$package_verification_code" \
  --arg supplier "Organization: The HDF Group" \
  --arg version "$version_info" \
  '{
    SPDXID: $spdxid,
    checksums: [
      {algorithm: "SHA256", checksumValue: $checksum}
    ],
    downloadLocation: $download_location,
    externalRefs: [
      {
        referenceCategory: "PACKAGE_MANAGER",
        referenceLocator: $purl,
        referenceType: "purl"
      }
    ],
    filesAnalyzed: true,
    licenseConcluded: $license,
    name: $name,
    originator: $originator,
    packageFileName: $package_file_name,
    packageVerificationCode: {
      packageVerificationCodeValue: $verification_code
    },
    primaryPackagePurpose: "SOURCE",
    supplier: $supplier,
    versionInfo: $version
  }' > "$package_json"

jq -cn \
  --arg relatedSpdxElement "$package_spdx_id" \
  '{
    relatedSpdxElement: $relatedSpdxElement,
    relationshipType: "DESCRIBES",
    spdxElementId: "SPDXRef-DOCUMENT"
  }' >> "$relationships_jsonl"

mkdir -p "$(dirname "$output_sbom")"

jq -n \
  --arg spdx_version "$spdx_version" \
  --arg document_spdx_id "SPDXRef-DOCUMENT" \
  --arg document_name "$document_name" \
  --arg document_namespace "$document_namespace" \
  --arg created "$created_at" \
  --arg creator_person "Person: $creator_person" \
  --arg creator_tool "Tool: $creator_tool" \
  --arg license_list_version "$license_list_version" \
  --arg data_license "$data_license" \
  --slurpfile packages "$package_json" \
  --slurpfile files "$files_jsonl" \
  --slurpfile relationships "$relationships_jsonl" \
  '{
    spdxVersion: $spdx_version,
    SPDXID: $document_spdx_id,
    name: $document_name,
    documentNamespace: $document_namespace,
    creationInfo: {
      created: $created,
      creators: [
        $creator_person,
        $creator_tool
      ],
      licenseListVersion: $license_list_version
    },
    dataLicense: $data_license,
    packages: $packages,
    files: $files,
    relationships: $relationships
  }' > "$output_sbom"

printf 'Wrote SPDX JSON SBOM to %s\n' "$output_sbom"
