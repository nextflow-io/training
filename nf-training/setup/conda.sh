# install conda
wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -f -p $HOME/miniconda
rm Miniconda3-latest-Linux-x86_64.sh
export PATH=$PATH:$HOME/miniconda/bin
echo "export PATH=$PATH:$HOME/miniconda/bin" >> ~/.bashrc
