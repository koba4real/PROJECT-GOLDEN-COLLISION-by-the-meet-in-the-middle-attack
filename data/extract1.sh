# Check if the correct number of arguments is passed
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <input_file.rtf>"
  exit 1
fi

# Input RTF file
input_file="$1"

# Temporary file for processing
temp_file="$(mktemp)"

# Extract and save the required lines using grep and regex
cat "$input_file" |
  grep -E "(Running with n=[^,]+,|Dictionary size: [^,]+|Root process execution time: [^,]+|Fill:[^,]+|Probe:[^,]+|Solution(s)? found:.*)" > "$temp_file"

# Overwrite the original file with the filtered output
mv "$temp_file" "$input_file"

echo "Filtered output saved to $input_file"