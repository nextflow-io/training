#!/usr/bin/env bash

# Customise the terminal command prompt
printf "export PS1='\\[\\e[3;36m\\]\${PWD#/workspaces/} ->\\[\\e[0m\\] '\n" >> $HOME/.bashrc
export PS1='\[\e[3;36m\]${PWD#/workspaces/} ->\[\e[0m\] '

# Update Nextflow
nextflow self-update
nextflow -version

# Build nf-test from main branch with v2 parser support (PR #336)
# TODO: Remove this once nf-test releases a version with strict syntax support
# See: https://github.com/askimed/nf-test/pull/336
NFTEST_COMMIT="350bb147a23a7f0aa657c13342d9726c0e3edacc"
echo "Building nf-test from commit ${NFTEST_COMMIT} (v2 parser support)..."

# Install Maven if not present
if ! command -v mvn &> /dev/null; then
    echo "Installing Maven..."
    apt-get update -qq && apt-get install -y -qq maven > /dev/null
fi

cd /tmp
git clone --quiet https://github.com/askimed/nf-test.git
cd nf-test
git checkout --quiet ${NFTEST_COMMIT}
mvn install -DskipTests -q
cp target/nf-test.jar ~/.nf-test/nf-test.jar
cd /workspaces/training
rm -rf /tmp/nf-test ~/.m2/repository
echo "nf-test updated with v2 parser support"

cat /usr/local/etc/vscode-dev-containers/first-run-notice.txt
