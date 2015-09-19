#!/bin/sh
printcore -f conf_power_map -v coefs_030 -c 0.5 1.5 0 -shift 4 0 -dbb 0 0 50 0
ps2pdf mycore1.eps _mycore1.pdf
pdfcrop _mycore1.pdf mycore1.pdf
rm _mycore1.pdf
