#!/usr/bin/env bash
set -euo pipefail
# test-split.sh
# Copyright Jackson M. Tsuji, 2025
# End-to-end test for ampliwrangler split

# GLOBAL variables
readonly SCRIPT_NAME="${0##*/}"
FAILED_TESTS=0

#### FUNCTIONS
#######################################
# Compares two MD5 hashes for the same file in HARD-CODED output dirs
# Globals:
#   FAILED_TESTS: count of failed tests
# Arguments:
#   input_basename: the basename (no directory, no extension) of the desired MD5 hashes to test between the outputs and expected_outputs folders
# Returns:
#   updates FAILED_TESTS
#######################################
function check_md5_files() {
  # User-provided inputs
  local test_filepath
  test_filepath=$1
  local reference_filepath
  reference_filepath=$2

  md5_test=$(cut -d " " -f 1 "${test_filepath}")
  md5_reference=$(cut -d " " -f 1 "${reference_filepath}")

  if [[ "${md5_test}" != "${md5_reference}" ]]; then
    echo "[ $(date -u) ]: '${test_filepath}': FAILED: Actual hash does not match expected hash"
    # TODO - change to local variable
    FAILED_TESTS=$((${FAILED_TESTS}+1))
  else
    echo "[ $(date -u) ]: '${test_filepath}.tsv': PASSED"
  fi
}
##########

# If no input is provided, provide help and exit
if [[ $# -eq 0 ]]; then
  echo "No arguments provided. Please run '-h' or '--help' to see help. Exiting..." >&2
  exit 1
elif [[ $1 = "-h" ]] || [[ $1 = "--help" ]]; then

  # Help statement
  printf "${SCRIPT_NAME}: run end-to-end test for ampliwrangler split\n"
  printf "Usage: ${SCRIPT_NAME} test_dir\n\n"
  printf "Positional arguments (required):\n"
  printf "   test_dir: path to the test directory containing the 'inputs' and 'outputs-expected' test folders\n\n"
  printf "Note: script will give an exit status of 1 if any tests fail; otherwise exit status will be 0.\n\n"

  # Exit
  exit 0
fi

### Get user inputs
test_dir=$1

echo "[ $(date -u) ]: Running ${SCRIPT_NAME}"
echo "[ $(date -u) ]: Test dir: '${test_dir}'"

### Expected positions of folders
input_dir="${test_dir}/inputs"
expected_outputs_dir="${test_dir}/outputs-expected"
output_dir="${test_dir}/outputs" # should NOT yet exist!

### Look for dirs
if [[ ! -d "${input_dir}" ]]; then
  echo "[ $(date -u) ]: ERROR: did not find input dir at '${input_dir}'. Exiting..."
  exit 1
elif [[ ! -d "${expected_outputs_dir}" ]]; then
  echo "[ $(date -u) ]: ERROR: did not find expected outputs dir at '${expected_outputs_dir}'. Exiting..."
  exit 1
fi

### Make sure the output dir doesn't exist
if [[ -d "${output_dir}" ]]; then
  echo "[ $(date -u) ]: ERROR: output directory '${output_dir}' already exists. Please specify an output dir that does not exist. Exiting..."
  exit 1
fi

## Test 1
test_ID="01_standard"
echo "[ $(date -u) ]: Running script on test data '${test_ID}'"
ampliwrangler split \
  -i "${input_dir}/manifest.tsv" \
  -o "${output_dir}" \
  -r "run-id" \
  -v \
  > "${test_dir}/${test_ID}.log" 2>&1
mv "${test_dir}/${test_ID}.log" "${output_dir}"
# TODO - print an error message if this fails

# Generate MD5 hash on expected output and compare
output_tsv_filepaths=($(find "${output_dir}" -name "*.tsv"))
for output_tsv_filepath in "${output_tsv_filepaths[@]}"; do
  output_tsv_filename="${output_tsv_filepath##*/}"
  md5sum "${output_tsv_filepath}" > "${output_tsv_filepath}.md5"
  check_md5_files "${output_tsv_filepath}.md5" "${expected_outputs_dir}/${output_tsv_filename}.md5"
done

### Overall status
if [[ "${FAILED_TESTS}" -eq 0 ]]; then
  echo "[ $(date -u) ]: All tests PASSED. Deleting output folder."
  rm -r "${output_dir}"
  echo "[ $(date -u) ]: Testing complete."
  exit 0
else
  echo "[ $(date -u) ]: ${FAILED_TESTS} test(s) FAILED. See above log for details. Keeping outputs dir for reference. Testing complete."
  exit 1
fi
