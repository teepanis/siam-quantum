# Compiling Guide

Siam Quantum executable can be downloaded for <span style="color: blue;">**Windows**</span>, <span style="color: blue;">**Linux**</span>, and <span style="color: blue;">**MacOS**</span> platform. However, you may want to compile it yourself for your code development project.

There are 2 dependencies: **LAPACK** and **LibXC** libraries. This guide is a __walkthrough__ on how to compile these libraries and Siam Quantum 2.

<br>

## <span style="color: blue;">Linux Platform</span>

We tested the walkthrough on Ubuntu 24.04.3 LTS AMD 64-bit distribution, which actually ran on a VirtualBox machine. So we included here exactly how we installed the Ubuntu on the VirtualBox. You may skip step 1) if you already have your Linux distribution working.

### <span style="color: green;">1) Install Ubuntu (on VirtualBox)</span>

- On a Windows 10 machine, we downloaded the Ubuntu ISO from its website, and saved the file: `ubuntu-24.04.3-desktop-amd64.iso`

- We opened VirtualBox application, and click `New` to setup a new virtual machine

- We specified the following information
> - VM Name: Ubuntu 24.04.3 LTS AMD 64-bit Desktop
> - VM Folder: C:\Users\Admin\VirtualBox VMs
> - ISO Image: >> just downloaded
> - OS: Linux
> - OS Distribution: Ubuntu
> - Proceed with Unattended Installation
> - RAM 4096, CPU 2, Disk Size 25GB
> - Username: vbox
> - Password: F=ma

- Click `Finish`, and the VirtualBox will start Ubuntu installation process. Be sure to click _Install Ubuntu_ afterward. When it finished copying all the files, it restarted the virtual machine automatically.

- We now have a freshly installed version of Ubuntu 24.04.3 (running on a virtual machine).

### <span style="color: green;">2) Prepare Ubuntu</span>

- We started up the Ubuntu, and prepared the development tools by executing the following commands.

```
$ sudo apt instal vim
$ sudo apt install build-essential
$ sudo apt install gfortran
$ sudo apt-get -y install cmake
$ gcc --version
```

- The last command showed us that the gcc version is `13.3.0`.

### <span style="color: green;">3) Compile LAPACK</span>

- We downloaded lapack-3.12.1.tar.gz from its website using Mozilla (inside Ubuntu) and saved it to `/home/vbox/Downloads` directory.

- We executed the following commands:

```
$ cd /home/vbox
$ mkdir local
$ cd /home/vbox/local
$ gzip -cd /home/vbox/Downloads/lapack-3.12.1.tar.gz | tar -xvf -
```

- The commands above should unpack the LAPACK nicely on your `$HOME/local`, we will then proceed to compile it. 

```
$ cd lapack-3.12.1
$ cp INSTALL/make.inc.gfortran make.inc
$ make
```

- We now have the library `liblapack.a` and `labrefblas.a` needed to compile Siam Quantum.


### <span style="color: green;">4) Compile LibXC</span>

As of 2026, the LibXC already released version 7, but we only tested Siam Quantum with version 5, which is 5.2.3 to be exact.

- We downloaded libxc-5.2.3.tar.gz and unpack it onto the directory `$HOME/local`.

- Here are the commands we used to compile it.

```
$ cd libxc-5.2.3
$ cmake -DCMAKE_INSTALL_PREFIX:PATH=/home/vbox/local/libxc-libs/ -H. -Bobjdir
$ cd objdir
$ make clean
$ make
$ make install
```
- Take note of the files in the directory `$HOME/local/libxc-libs`. We will need them for compiling Siam Quantum.
- You may test if the LibXC has been compiled successfully by executing the program `xc-info` in `$HOME/local/libxc-libs/bin`. If it prints usage information, then everything went well so far.

### <span style="color: green;">5) Compile Siam Quantum</span>

- Download the source code into `$HOME/local`, unpack it, and inspect its directory structure.
```shell
$ cd $HOME/local
$ unzip sq-2.0.0.zip
$ cd sq-2.0.0
$ ls
CHANGELOG.md   README.md    basis    docs    examples    src
$ cd src
$ ls
Makefile basis.c basis.h check.c check.h ... libs
```
- From this output you see there is a directory: `sq-2.0.0/src/libs`. This is where we will copy the LAPACK and LibXC libraries into, using these commands:
```shell
$ cd libs

$ cp $HOME/local/lapack-3.12.1/liblapack.a ./
$ cp $HOME/local/lapack-3.12.1/liblrefblas.a ./

$ cp $HOME/local/libxc-libs/lib/libxc.a ./
$ cp $HOME/local/libxc-libs/include/*.h ./
```
- You should now see these files in the directory `sq-2.0.0/src/libs`:
```shell
ls $HOME/local/sq-2.0.0/src/libs
liblapack.a    libxc.a   xc_funcs.h           xc_funcs_worker.h
librefblas.a   xc.h      xc_funcs_removed.h   xc_version.h
```
- Edit the `src/Makefile`, and __un-comment__ the correct compiler version or edit it to reflect your gcc version. For example,
```Makefile
####### BEGIN SELECT COMPILER #######

## Windows MinGW using gcc compiler
#CC  = gcc
#CCOPT = -O3 -Wall -DLIBXC -I${INC_PATH} 
#LDOPT = -O3 -Wall

## GNU Linux using gcc compiler
CC  = gcc
CCOPT = -O3 -m64 -Wall -msse2 -DLIBXC -I${INC_PATH}
LDOPT = -O3 -m64 -Wall

## macOS using homebrew compiler
#CC = /opt/homebrew/bin/gcc-15
#CCOPT = -O3 -m64 -Wall -DLIBXC -I${INC_PATH}
#LDOPT = -O3 -m64 -Wall

######### END SELECT COMPILER #######
```
- Finally, we build Siam Quantum using these commands:

```shell
cd $HOME/sq-2.0.0/src
make clean
make
```
- We should now have the executable `sq`. Congratulation 🎉😃.


<br>

## <span style="color: blue;">Windows Platform</span>

We tested the the compilation on Intel(R) Core(TM) i7-9700 CPU @ 3.00GHz 3.00 GHz with Windows 10. The compilation platform for this test is `MinGW`.

### <span style="color: green;">1) Install MinGW</span>
- Download installation file from its SourceForge website: https://sourceforge.net/projects/mingw/files/Installer/mingw-get-setup.exe
- Double click it, and choose `Install`. Select the default directory
`C:\MinGW`, continue until this initial process is complete.
- In `MinGW Installation Manager` we are to select packages to install. Focus on the _Basic Setup_ tab, and select all of these:
> - mingw-developer-toolkit,
> - mingw32-base
> - mingw32-gcc-ada
> - mingw32-gcc-fortran
> - mingw32-gcc-g++
> - mingw32-gcc-obj
> - mysys-base  
- Now focus on _All Packages_ tab, search for these entries, and select them.
> - msys-unzip bin
> - msys-unzip doc
> - msys-unzip lic
- On the upper menu, choose Installation >> Apply Changes
- Open a terminal for the MinGW environment by browsing your file system through `C:\MinGW\msys\1.0`. You will find `msys.bat` here. Double click it.
- A shell terminal with the prompt `$` appears. To understand its file system we can execute:

```shell
$ pwd
$ mkdir local
$ gcc --version
```

- From above, I saw that the current directory is `/home/Admin`. Yours might not be `Admin` but something else, depending on the Windows username you are using. MinGW also reported that my gcc version was 6.3.0.
- _Under Windows_, use your mouse to browse your files through `C:\MinGW\msys\1.0\home\Admin`, you will see a directory `local` you have just created here.
- My point is this. A file in `/home/Admin` under _MinGW_ is the same file in `C:\MinGW\msys\1.0\home\Admin` _under Windows_. Hence, you can easily move a file in and out between the MinGW platform and the Windows platform.

### <span style="color: green;">2) Compile LAPACK</span>

- Use your favorite Windows web-browser like Explorer or Chrome to download lapack-3.12.1.tar.gz then move the file from your Downloads directory to `/home/Admin/local` under the MinGW platform. For example, use your mouse to drag the file over to the `C:\MinGW\msys\1.0\home\Admin\local` directory.
- Open the MinGW terminal (by double-click `C:\MinGW\msys\1.0\msys.bat`) and execute the following commands to unpack and compile LAPACK.


```shell
$ cd local
$ gzip -cd lapack-3.12.1.tar.gz | tar -xvf -
$ cd lapack-3.12.1
$ cp INSTALL/make.inc.gfortran make.inc
$ make
```

- We now have the library `liblapack.a` and `librefblas.a` needed to compile Siam Quantum. _Notes: The compilation failed at python TESTING, but the libraries were created just fine._

### <span style="color: green;">3) Compile LibXC</span>

- Again, using Windows web-browser to download libxc-5.2.3.tar.gz, then move it from your Downloads directory to `/home/Admin/local` under the MinGW platform.
- Open the MinGW terminal and proceed with the following commands. (Be sure to change `Admin` to your Windows username)

```shell
$ cd local
$ gzip -cd libxc-5.2.3.tar.gz | tar -xvf -
$ cd libxc-5.2.3
$ autoreconf -i
$ ./configure --prefix=/home/Admin/local/libxc-libs
$ make
$ make install
```

- Watch out for the _occasional_ error after unpacking the file (i.e. gzip -cd ...). Should that happens, repeat the unpacking again until all the files are unpacked without error, before proceeding to compile LibXC.
- At this point, we should have the LibXC library located in `/home/Admin/local/libxc-libs`. To test if the compilation is successful, execute:

```
$ /home/Admin/local/libxc-libs/bin/xc-info.exe
```

- If you see the MinGW terminal printing usage information, then you are on the right track. You can even use your mouse to drag the `xc-info.exe` to your Desktop folder and use Windows Command Prompt to execute it. If successful, as it should, it confirms that the executable file is Windows native, meaning it can be executed within native Windows environment.

### <span style="color: green;">4) Compile Siam Quantum</span>

- Unpack the source code distribution onto `$HOME/local`, and inspect its directory structure.

```shell
$ cd $HOME/local
$ unzip sq-2.0.0.zip
$ cd sq-2.0.0
$ ls
CHANGELOG.md   README.md    basis    docs    examples    src
$ cd src
$ ls
Makefile basis.c basis.h check.c check.h ... libs
```

- From this output you see there is a directory: `sq-2.0.0/src/libs`. This is where we will copy the LAPACK and LibXC libraries into, using these commands:

```shell
$ cd libs

$ cp $HOME/local/lapack-3.12.1/liblapack.a ./
$ cp $HOME/local/lapack-3.12.1/liblrefblas.a ./

$ cp $HOME/local/libxc-libs/lib/libxc.a ./
$ cp $HOME/local/libxc-libs/include/*.h ./
```

- You should now see these files in the directory `sq-2.0.0/src/libs`:

```shell
ls $HOME/local/sq-2.0.0/src/libs
liblapack.a    libxc.a   xc_funcs.h           xc_funcs_worker.h
librefblas.a   xc.h      xc_funcs_removed.h   xc_version.h
```

- Edit the `src/Makefile`, and __un-comment__ the correct compiler version or edit it to reflect your homebrew compiler version. For example,

```Makefile
####### BEGIN SELECT COMPILER #######

## Windows MinGW using gcc compiler
CC  = gcc
CCOPT = -O3 -Wall -DLIBXC -I${INC_PATH} 
LDOPT = -O3 -Wall

## GNU Linux using gcc compiler
#CC  = gcc
#CCOPT = -O3 -m64 -Wall -msse2 -DLIBXC -I${INC_PATH}
#LDOPT = -O3 -m64 -Wall

## macOS using homebrew compiler
#CC = /opt/homebrew/bin/gcc-15
#CCOPT = -O3 -m64 -Wall -DLIBXC -I${INC_PATH}
#LDOPT = -O3 -m64 -Wall

######### END SELECT COMPILER #######
```

- Finally, we build Siam Quantum using these commands:

```shell
cd $HOME/sq-2.0.0/src
make clean
make
```

- We should now have the executable `sq.exe`. Congratulation 🎉😃.

<br>

## <span style="color: blue;">MacOS Platform</span>

We tested the compilation on Macbook Pro M1 using _homebrew_ as the compiling platform.  

### <span style="color: green;">1) Prepare Homebrew</span>

- If not already, you need to install homebrew platform on your macos. Please consult its website: https://brew.sh/ to install it. 

- Open Terminal and execute the following commands to check if the homebrew packages are ready.

```shell
$ brew list
$ brew install vim
$ brew install cmake
$ brew install gcc
$ brew list gcc
```

- We will be using _vim_, _cmake_ and _gcc_ packages; so the command above is to make sure we have all the packages we need. In particular, the last command tells us where the `gcc` is located. Take note of the output. For example, for our case, we have the following output:

```
/opt/homebrew/Cellar/gcc/15.2.0/bin/gcc-15
/opt/homebrew/Cellar/gcc/15.2.0/bin/gfortran-15
```

### <span style="color: green;">2) Compile LAPACK</span>

- Download and save lapack-3.12.1.tar.gz in your `$HOME/local`, unpack it, and prepare the Makefile using the following commands:

```shell
$ cd $HOME/local
$ gzip -cd lapack-3.12.1.tar.gz | tar -xvf -
$ cd lapack-3.12.1
$ cp INSTALL/make.inc.gfortran make.inc
$ vim make.inc
```

- In the last command, we use _vim_ text editor to edit the file, specifying which compilers to use. Or use other text editors you feel comfortable with. In the file `make.inc`, look for the following entries and change their values according to the homebrew compilers listed above. For example,

```
CC = /opt/homebrew/Cellar/gcc/15.2.0/bin/gcc-15
FC = /opt/homebrew/Cellar/gcc/15.2.0/bin/gfortran-15
```

- Now we are ready to compile LAPACK, simply execute:

```
make
```

- We now have the library `liblapack.a` and `librefblas.a` needed to compile Siam Quantum.


### <span style="color: green;">3) Compile LibXC</span>

- Using web-browser like Chrome or Safari to download libxc-5.2.3.tar.gz, then move it from your Downloads directory to `$HOME/local`. If you are a new mac user, try these commands:

```
$ ls $HOME/Downloads/
$ cp $HOME/libxc-5.2.3.tar.gz  $HOME/local
```

- Proceed to unpack it, and use cmake to compile it.

```shell
$ cd libxc-5.2.3
$ cmake -DCMAKE_INSTALL_PREFIX:PATH=$HOME/local/libxc-libs/ -H. -Bobjdir  \
-DCMAKE_C_COMPILER=/opt/homebrew/Cellar/gcc/15.2.0/bin/gcc-15  \
-DCMAKE_Fortran_COMPILER=/opt/homebrew/Cellar/gcc/15.2.0/bin/gfortran-15  \
-DCMAKE_POLICY_VERSION_MINIMUM=3.5
$ cd objdir
$ make clean
$ make
$ make install
```

- The 2nd command is quite long. It tells _cmake_ which C and Fortran compiler to use; and the last flag, DCMAKE_POLICY_VERSION_MINIMUM, is to override LibXC's cmake requirement. 

- At this point, we should have the LibXC library located in `$HOME/local/libxc-libs`. To test if the compilation is successful, execute:

```shell
$ ./xc-info
```

- If you see usage information, then you are on the right track.


### <span style="color: green;">4) Compile Siam Quantum</span>

- Unpack the source code distribution onto `$HOME/local`, and inspect its directory structure.

```shell
$ cd $HOME/local
$ unzip sq-2.0.0.zip
$ cd sq-2.0.0
$ ls
CHANGELOG.md   README.md    basis    docs    examples    src
$ cd src
$ ls
Makefile basis.c basis.h check.c check.h ... libs
```

- From this output you see there is a directory: `sq-2.0.0/src/libs`. This is where we will copy the LAPACK and LibXC libraries into, using these commands:

```shell
$ cd libs

$ cp $HOME/local/lapack-3.12.1/liblapack.a ./
$ cp $HOME/local/lapack-3.12.1/liblrefblas.a ./

$ cp $HOME/local/libxc-libs/lib/libxc.a ./
$ cp $HOME/local/libxc-libs/include/*.h ./
```

- You should now see these files in the directory `sq-2.0.0/src/libs`:

```shell
ls $HOME/local/sq-2.0.0/src/libs
liblapack.a    libxc.a   xc_funcs.h           xc_funcs_worker.h
librefblas.a   xc.h      xc_funcs_removed.h   xc_version.h
```

- Edit the `src/Makefile`, and __un-comment__ the correct compiler version or edit it to reflect your homebrew compiler version. For example,

```Makefile
####### BEGIN SELECT COMPILER #######

#CC  = icc
#CCOPT = -O3 -funroll-loops -fno-alias -m64 -Wall
#LDOPT = -O3 -m64 -Wall

## GNU Linux using gcc compiler
#CC  = gcc
#CCOPT = -O3 -m64 -Wall -msse2 -DLIBXC -I${INC_PATH}
#LDOPT = -O3 -m64 -Wall

## macOS using homebrew compiler
CC = /opt/homebrew/bin/gcc-15
CCOPT = -O3 -m64 -Wall -DLIBXC -I${INC_PATH}
LDOPT = -O3 -m64 -Wall

######### END SELECT COMPILER #######
```

- Finally, we build Siam Quantum using these commands:

```shell
cd $HOME/sq-2.0.0/src
make clean
make
```

- We should now have the executable `sq`. Congratulation 🎉😃.

