Overview
========

DoGPack stands for Discontinuous Galerkin Package.
It is a set of C++ libraries intended to be used by
applications to solve hyperbolic partial differential equations.
A typical dogpack application solves a hyperbolic equation of the
form 

       q_t + div(f) = \psi.

The user is expected to supply the flux
function f and the source term \psi as library callbacks.

Installation instructions are found at the following website

       http://www.dogpack-code.org/install.html

To obtain a fresh clone of the repository, use:

    git clone https://[username]@bitbucket.org/imsejae/dogpack-developer.git dogpack.git

  or for faster transfer (requires .ssh configuration described below)

    git clone git@bitbucket.org:imsejae/dogpack-developer.git

see `Installation Instructions' for further details.

Class Hierarchy
===============

DoGPack uses a class hierarchy for every 'application', which is a problem the
user wishes to run.  Similar to Clawpack [http://depts.washington.edu/clawpack/], 
all applications live under the apps tree, and common library routines live in lib.

The following tree describes the hierarchy of the 2D structured examples

  lib/DogSolver <- DogSolverTB (inside DogSolver.h) <- DogSolverCart2 <- AppSolverCart2 <- AppSolver,

where '<-' means inherits from.

DogSolver and DogSolverTB classes provide a handful of pure virtual methods which must be
defined by each derived class.  For example, one such function is,

  virtual double GetCFL(double dt)const=0;

Without the '=0' at the end, a derived class may choose to overide it, whereas
in this case, each derived class must define each of these virtual methods.

Alongside each of the 'Solver' classes are a hierarchy of 'State' classes which
contain information about all dynamic state variables.  The idea of casting
state variables into a class means that users should be able to revert back to
previous states of the system, for example for restarts or rejecting time
steps, and class encapsulation means that all relevant dynamic state variables
should be self contained.  For example, the following describes a 2D-structured
example,

lib/DogState <- DogStateTB (inside DogState.h) <- DogStateCart2 <- AppStateCart2

DogSolver Class
---------------

TODO - write some notes on what the purpose of DogSolver is.  Better yet, we
should be writing more notes inside the header files ... (-DS)

Installation Instructions
==================================

  To use the git transfer protocol you need to configure your ~/.ssh
  subdirectory properly as follows (see step 6 of
  "https://confluence.atlassian.com/pages/viewpage.action?pageId=270827678"):

    === need the following lines in your ~/.ssh/config file: ===
    Host bitbucket.org
     User [username]
     IdentityFile ~/.ssh/id_rsa

  and of course this means you need your
  ~/.ssh/id_rsa and ~/.ssh/id_rsa.pub files set up,
  as described in "man ssh" and "man ssh-keygen":

    ssh-keygen

  Then copy and paste the contents of your new
  ~/.ssh/id_rsa.pub into bitbucket's list of your SSH keys,
  accessible after logging into your bitbucket account
  by clicking on the tab with your name, selecting "Account",
  selecting "SSH keys" under "Account Settings", and pasting
  your public key into the SSH Key box.

The first option is to work with the repository directly using git.

  A good source of information about git is the Pro Git book:

    http://git-scm.com/book

  Before committing anything with git, you should edit your ~/.gitconfig file
  and ensure that your name and email are configured:

     [user]
     name  = [username]
     email = [yourname@email.com]

An alternative interface to git is supplied by the easy-git (eg) script:

    cd ~/bin # or wherever you want to put the shell script
    # download the latest version from http://people.gnome.org/~newren/eg/download/
    curl -o eg http://people.gnome.org/~newren/eg/download/1.7.3/eg
    chmod +x eg

  easy-git (eg) smooths the warts off the git interface and provides
  a svn-like interface.  eg download and documentation is available at:

    http://people.gnome.org/~newren/eg/

A third option is to use mercurial (hg) to work with the repository.

  A good source of information about hg is available at:

    http://hgbook.red-bean.com/

  If you wish to use hg to push to or pull from the git repository,
  you can do so by using hg-git.

  Information is available from

    http://hg-git.github.com

  and download is available from

    https://bitbucket.org/durin42/hg-git

  Detailed information on installing and using hg-git follows:

    # install hg-git using something like one of the following:
    easy_install hg-git # fails on my mac. -eaj
    fink update hg-git-py27 # worked on my mac. -eaj
    # edit ~/.hgrc to enable hg-git by adding two lines to [extensions] section:
      [extensions]
      hgext.bookmarks =
      hggit =
    # if we were using github we could simply do the following:
    #hg clone https://eajohnson@bitbucket.org/imsejae/dogpack-git.git dogpack.hg
    # and then we could push and pull directly from the local hg repository.
    # Unfortunately bitbucket does not support hg-git, so we must
    # create a bare git repository to serve as an intermediary:
    git clone --mirror https://eajohnson@bitbucket.org/imsejae/dogpack-developer.git dogpack.mirror
    # clone the bare repository
    hg clone dogpack.bare dogpack.hg

    # work with the hg clone:
    cd dogpack.hg
    [... do some work ...]
    hg commit
    # transfer the work to the intermediary
    hg push
    # transfer the work to the source tree
    cd ../dogpack.bare
    git push

    # update the hg clone from the repository
    git fetch --all
    cd ../dogpack.hg
    hg pull
    hg update

  By using hg hooks (e.g. the "outgoing" and "pre-pull"
  hooks you should be able to configure hg so that
  hg push and hg pull automatically cause "git push"
  and "get fetch --all" to be called in the local
  intermediary git repository.

