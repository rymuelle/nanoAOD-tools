#!/bin/bash
export X509_USER_PROXY="$PWD/x509up_u63115"
voms-proxy-init -voms cms
