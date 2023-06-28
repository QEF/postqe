# !/bin/env bash
#

curl -SL https://github.com/QEF/q-e/archive/refs/tags/qe-7.2.tar.gz -o qe-7.2.tar.gz
echo `md5sum qe-7.2.tar.gz` | grep --quiet "^ece7dcde72915110cafc7ae3c7c99da3[[:blank:]]"

