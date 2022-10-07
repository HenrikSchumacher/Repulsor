# Repulsor

# Installation

After cloning make sure to run the following to connect all submodules to their repos.

    git submodule update --init --recursive
    

Pull changes from the remote repositories of any submodule by executing

    git submodule update --remote --recursive
    
    
# Trouble shooting

If you accidentally modified one of the submodules you can run

    git submodule foreach --recursive git reset --hard
    
to repair this.

