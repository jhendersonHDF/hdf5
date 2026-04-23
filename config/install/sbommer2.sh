spdx_safe_id_suffix() {
  printf '%s' "$1" | sed -E 's/[^A-Za-z0-9.-]+/-/g; s/^-+//; s/-+$//'
}

# MAIN LOGIC

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

# Create 'package' JSON file
jq -cn ...

# Append to 'relationships' file
jq -cn \
  --arg relatedSpdxElement "$package_spdx_id" \
  '{
    relatedSpdxElement: $relatedSpdxElement,
    relationshipType: "DESCRIBES",
    spdxElementId: "SPDXRef-DOCUMENT"
  }' >> "$relationships_jsonl"

# Assemble final output file
jq -n ...
