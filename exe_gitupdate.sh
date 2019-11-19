#!/bin/bash
#$ -S /bin/bash
#$ -cwd

sudo git pull
#sudo chmod -R 775 *

sudo git add --all
sudo git commit -m 'Minor Update'
sudo git push origin master
