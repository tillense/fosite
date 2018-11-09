#!/usr/bin/env bash

echo "test", $1

awk '{gsub("<div class=\"textblock\">", "")}1' ${1}/doc/html/citelist.html > ${1}/doc/html/citelist_tmp.html
mv ${1}/doc/html/citelist_tmp.html ${1}/doc/html/citelist.html
sed -i -- 's/href="index.html"/href="main.html"/g' ${1}/doc/html/*.html
sed -i -- 's/{text:"Main",url:"index.html"}/{text:"Main",url:"main.html"}/g' ${1}/doc/html/menudata.js
mv ${1}/doc/html/index.html ${1}/doc/html/main.html
sed -i -- 's/href="fosite.html"/href="index.html"/g' ${1}/doc/html/*.html
mv ${1}/doc/html/fosite.html ${1}/doc/html/index.html
sed -i -- 's/class="code"//g' ${1}/doc/html/*.html
