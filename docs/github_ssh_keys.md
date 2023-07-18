Set up SSH key for GitHub
=========================

For more advice: https://docs.github.com/en/authentication/connecting-to-github-with-ssh/about-ssh

After checking for existing keys, if you receive error that ~/.ssh doesn't exist then you don't have one. If there already is one (ie. id_rsa.pub, id_ed25519.pub) then you can either connect it to GitHub or generate new one.
```
ls -al ~/.ssh #check for existing keys
ssh-keygen -t ed25519 -C "your_email@example.com"                                       #use your GitHub email address
#Enter a file in which to save the key (/c/Users/you/.ssh/id_algorithm):[Press enter]
#Enter passphrase (empty for no passphrase): [Type a passphrase]
eval "$(ssh-agent -s)"                                                                  #start ssh-agent
ssh-add ~/.ssh/id_ed25519                                                               #add your SSH private key to ssh-agent
clip < ~/.ssh/id_ed25519.pub                                                            #copy SSH public key 
```
After copying your SSH public key, go to GitHub --> Settings --> SSH and GPG keys (under Access) --> Add new public SSH key

To test connection 
```
ssh -T git@github.com 
```
A successful connection should result in 
> Hi username! You've successfully authenticated, but GitHub does not provide shell access.

Activate the environment
```