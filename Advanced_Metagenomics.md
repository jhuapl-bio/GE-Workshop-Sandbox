# Advanced Metagenomics Installation

## WSL2 (Windows Only)

Prerequisites: None

## Nextflow 

Prerequisites: 

1. Java (OpenJDK)
2. Unix Environment (WSL2 see above for Windows Users)


On Windows, you need to be sure you are running WSL2

## Docker

Prerequisites: None (But you need Admin Permissions)

[Mac OSX](https://docs.docker.com/desktop/install/mac-install/)

Make sure to select the right processor type. Typically, newer Mac Models are the Apple Silicon option

[Windows](https://docs.docker.com/desktop/install/windows-install/)

You will need admin access to install this

[Linux](https://docs.docker.com/engine/install/ubuntu/)

Make sure to follow the steps in the `post-installation` steps

Or you can follow these commands below. Simply copy and paste them into your terminal


```
# Add Docker's official GPG key:
sudo apt-get update
sudo apt-get install ca-certificates curl gnupg
sudo install -m 0755 -d /etc/apt/keyrings
curl -fsSL https://download.docker.com/linux/debian/gpg | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg
sudo chmod a+r /etc/apt/keyrings/docker.gpg

# Add the repository to Apt sources:
echo \
  "deb [arch="$(dpkg --print-architecture)" signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/debian \
  "$(. /etc/os-release && echo "$VERSION_CODENAME")" stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update


sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin


sudo groupadd docker
sudo usermod -aG docker $USER

## only run these commands once (between the ##) as repeat lines can mess up your env. Use nano if you need to edit it and this was run already
sudo sed -i "1s/^/$USER:$(id -u):1\n/" /etc/subuid
sudo sed -i "1s/^/$USER:$(id -g):1\n/" /etc/subgid
sudo apt install jq

echo $(jq --arg user "$USER" '. += {"userns-remap": $user}' /etc/docker/daemon.json) > ~/daemon.json && sudo mv ~/daemon.json /etc/docker/daemon.json

##

```


### Base Install

### Images

## Basestack

Prerequisites: Docker

1. Go to [here](https://github.com/jhuapl-bio/Basestack/releases/latest)
2. Download the latest binary 
    a. Mac OSX: `.dmg`
    b. Windows (non-admin): `win-unpacked.zip`. You will need to extract/unzip and double click the `.exe` each time to run
    c. Windows (admin): `*Setup.exe`
    d. Linux: `AppImage`. Make sure to select `x86_64` (most cases for your laptop)
        - You will need to run `chmod +x` on the AppImage to allow it to be double-clickable. Otherwise run with `./Basestack.x86_64.AppImage`

## Conda 

Prerequisites: None

Ensure that you are INSIDE `WSL2` (see steps above to install). When running the script, select `Enter` or type `yes` when prompted

1. `WSL2` or `Linux`: 

```
cd $HOME
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

2. `Mac OSX`

`ARM64`

```
cd $HOME
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
bash Miniconda3-latest-MacOSX-arm64.sh
```

`AMD64`

```
cd $HOME
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh
```

## Git


