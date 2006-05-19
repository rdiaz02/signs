#! /bin/bash
cd /http/signs/Rpackage.and.Data
R CMD build SurvSignature.01
sudo R CMD INSTALL SurvSignature.01_0.1.tar.gz
for machine in $(seq 1 30); do scp SurvSignature.01_0.1.tar.gz 192.168.2.$machine:~/.; done
for machine in $(seq 1 30); do ssh 192.168.2.$machine 'sudo R CMD INSTALL SurvSignature.01_0.1.tar.gz'; done

