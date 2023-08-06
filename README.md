# Repulsor

# Installation

Either clone with

    git clone --recurse-submodules

or clone as usual and then run the following to connect all submodules to their repos.

    git submodule update --init --recursive
    

Pull changes from the remote repositories of any submodule by executing

    git submodule foreach --recursive git checkout main
    git submodule foreach --recursive git pull
    
    
# Trouble shooting

If you accidentally modified one of the submodules you can run

    git submodule foreach --recursive git reset --hard
    
to repair this.

