---
layout: default
title: Introduction_to_Linux
rank: 1
---

# Introduction to Linux
### Author: Srikeerthana Kuchi
#### Contact: Srikeerthana.Kuchi@glasgow.ac.uk

```
Command line shortcuts 

Up/Down arrows: Previous commands
!!: Reruns previous command
Tab: Auto complete
Tab+Tab: All available options 
Ctrl+a: Move cursor to start of line
Ctrl+e: Move cursor to end of line
Alt+: Alternates between terminals
Ctrl+l: Clear screen (or Command+k on Mac) 
Ctrl+c: Terminates the running program 
Ctrl+z: Suspends the running program
Ctrl+w: Removes a previous word
Ctrl+d: Logout
Ctrl+d(in a command): Removes a character 
Ctrl+u: Removes till the beginning 
```

## Introduction

[Linux](#Info) is an open-source operating system [OS](#Info) developed based on the kernel created by Linus Benedict Torvalds. In the last two decades, Linux has gained much popularity and now is being used on many platforms. Nowadays, most high-end servers to mobile phones (Android OS or iOS) run on different variants of Linux.

Linux computers/servers are installed for multi-user usage. In this course, we will work on a high performance cluster machine running [Ubuntu](#Info) server edition. Most of the commands specified in this manual can be used in any other distribution (i.e., CentOS, Debian, etc.) of Linux operating system. 

#### Resource: How to install Ubuntu?

To install a desktop edition of Ubuntu on personal computers, please follow the instructions in the following link.
[Install Ubuntu](https://tutorials.ubuntu.com/tutorial/tutorial-install-ubuntu-desktop#0)

### The Terminal

We use terminal (AKA command line interface) to interact with the operating system. The terminal by default runs one of the “shells”. Shell is a program that sits between the user and the [kernel](#Info) and translates user commands (text) into machine code. The advantages of using command line are greater control and flexibility over the system or software and multiple commands can be saved in a file and executed as a program.

The most common shells are:

Bourne Shell 
Bourne Again Shell – BASH (variant is Z Shell) 
C Shell (variant is T Shell) 
K Shell 

Among these Bourne Again Shell (BASH) is the most popular one. This is the default shell on the system, and we will be using it throughout this course. 

#### Connecting to Linux Server

In this course, we will be using MobaXterm application to access the Ubuntu OS. Please use the provided username and password to login into your account.

> Open MobaXterm -> Sessions -> New session -> [ssh](#Info) -> add remote host and username -> OK -> Enter password -> Don't save password

> Files can be downloaded and uploaded into the server

### Linux command structure

When you open a terminal in Linux, (MobaXterm by default opens a terminal), you will see a command prompt, ready to take commands. The default location on the terminal is your “home directory”. It is represented with ~ (tilde) symbol.

Copy the command below and paste it into your command line to copy the contents of the directory Linux to your home directory.

```bash
cp -R /home4/VBG_data/Linux .
```

All Linux commands are single words (can be alpha-numeric), with optional parameters followed by arguments. For historical reasons, some of the early commands are only two letter long and case sensitive. Most of the command options (also called flags) are single letters. They should be specified after the command before giving any input. 

```bash
ls -l Linux
```

"ls" is the command to list the contents of the directory, "-l" is the option for long listing and "Linux" is the input, which is optional in this case. Without the input, "ls" shows all the contents of the current directory (Type <mark style="background-color: lightblue">ls -l</mark>). 

To clear the terminal screen,

```bash
clear
```

### First Commands

Directories are the Unix equivalent of folders on a [PC](#Info) or a [Mac](#Info). They are organized in a hierarchy, so directories can have sub-directories and so on. Directories, like folders, are useful to keep your data files organized. The location or directory you are currently in, is called the current working directory. The location or “full pathname” of the file SARS-CoV-2.fa in the 'Linux' directory can be expressed as: 

```
Do not type this - won't work
/home_location/username/Linux
```

#### Tab completion

Typing out longer file names can be boring, and you are likely to make typos that will, at best, make your command fail with a strange error and at worst, overwrite some of your carefully crafted analysis. 

Tab completion is a trick that normally reduces this risk significantly. Instead of typing out "ls Interesting_stuff/", try typing “ls Int” and press the Tab button (instead of Enter). The rest of the folder/file names that begin with “Int” should be listed. If you have two folders/files with similar names (e.g., my_awesome_scripts/ and my_awesome_results/) then you might need to give your terminal a bit of a hand to work out which one you want. In this case if you type "ls –l m", when you press Tab the terminal would read "ls –l my_awesome_". You could then type "s" followed by another press of Tab button and it would figure out that you meant "my_awesome_scripts/".

### Info

| Terminology | Description                                                                                 |
| ----------- | ------------------------------------------------------------------------------------------- |
| Linux       | Unix derivative, most popular variant of Unix                                               |
| OS          | Software that commands the hardware and make the computer work                              |
| Ubuntu      | Free Linux distribution (distro) based on Debian (an oldest OS based on Linux kernel)       |
| Kernel      | Core interface between a computer’s hardware and its processes, manages available resources |
| ssh         | Program for logging in to a remote machine specified with a host name                       |
| PC          | A personal computer                                                                         |
| Mac         | A Macintosh computer                                                                        |

***
### **Points to remember**:                                         
> Linux commands are case sensitive and are always single words \
> Options follow the command - and they start with a single hyphen (-) and a character or a double hyphen (- -) and a word \
> Single character options can be combined \
> Argument can be one or two inputs \
> You can write more than one command separating with a semicolon; You can use “tab” to auto-fill the command
***

### Important Commands

*(a)* ls \
Lists information about the files/directories. Default is the current directory. Sorts entries alphabetically. 

Commonly used options: 
-l long list 
-a show all files (including hidden files) 
-t sort based on last modified time 

```bash
ls -l
```

Information (from left to right): 
•    File permissions 
•    Number of links 
•    Owner name 
•    Group name 
•    Number of bytes 
•    Abbreviated month, last modified date and time 
•    File/Directory name 

*(b)* pwd \
Returns the path of the current working directory (print working directory) to the standard output. 

```bash
pwd
```

*(c)* cd \
Change current working directory to the specified directory. 

```bash
cd Exercises/
pwd
```

We are now in the directory "Exercises". Typing the command "cd .." changes it to the parent directory from which the previous command was typed in. Typing "cd" will change the current directory to the home directory.

```bash
cd
cd Linux/
```

*(d)* mkdir \
This command creates a directory in the current working directory if no file/directory exists with the specified name. 

```bash
mkdir Practice
ls -l
```

*(e)* rmdir \
This command is used to remove directories. 

```bash
rmdir Practice
ls -l
```

*(f)* touch \
It is file’s time-stamp changing command. However, it can be used to create an empty file. This command is generally used to check if there is write permission for the current user.

```bash
touch temp-file
ls -l
```

*(g)* rm \
rm is used for removing files and directories.

```bash
rm temp-file
ls -l
```

> [!WARNING]
> **To remove directories use "-r" option. Please remember once a file or directory is deleted, it will not go to "Recycle bin" in Linux and there is no way you can recover it.**

*(h)* cp \
 Copies the content of the source file/directory to the target file/directory. To copy directories, use "-r" option.

```bash
touch temp1
cp temp1 temp2
ls -l
```

*(i)* mv \
To move/rename a file or a directory.

```bash
mkdir temp
mv temp1 temp/.
mv temp2 temp3
ls -l
```

The second command moves the "temp1" file into the directory "temp". The "." (dot) at the end of the command retains the name of the file, whereas the third command renames the file "temp2" to "temp3".

*(j)* ln \
Link command is used to make links to files/directories. We encourage you to create links rather than copying data in order to save space.

```bash
ln -s temp/temp1 .
ls -l 
```

### File viewers

*(a)* cat \
The concatenate command combines files (sequentially) and prints on the screen (standard output).

```bash
cat SARS-CoV-2.fa
```

*(b)* more/less \
These commands are used for viewing the content of the files; faster with large input files than text editors; not the entire file is read at the beginning.

```bash
more SARS-CoV-2.fa
```

Press “Enter” to view lines further and “q” to quit the program

*(c)* head/tail \
These commands show first/last 10 lines (default) respectively from a file.

```bash
head SARS-CoV-2.fa
```

### File editors

There are many non-graphical text editors like ed, emacs, vi and nano available on most Linux distributions. Some of them are very sophisticated (e.g., vi) and for advanced users. 

Nano (earlier called pico) is like any graphical editor without a mouse. All commands are executed using the keyboard, using the <CTRL> key modifier. It can be used to edit virtually any kind of text file from the command line. Nano without a filename gives you a standard (blank) nano window. 

At the bottom of the screen, there are commands with a symbol in front. The symbol tells that you need to hold down the Control (Ctrl) key, and then press the corresponding letter of the command you wish to use. 

Ctrl+X will exit nano and return you to the command line. 

```
Nano quick reference

Ctrl+X: Exit the editor. If you’ve edited text without saving, you’ll be prompted as to whether you really want to exit. 

Ctrl+O: Write (output) the current contents of the text buffer to a file. A filename prompt will appear; press Ctrl+T to open the file navigator shown above. 

Ctrl+R: Read a text file into the current editing session. At the filename prompt, hit Ctrl+T: for the file navigator. 

Ctrl+K: Cut a line into the clipboard. You can press this repeatedly to copy multiple lines, which are then stored as one chunk. 

Ctrl+J: Justify (fill out) a paragraph of text. By default, this reflows text to match the width of the editing window. 

Ctrl+U: Uncut text, or rather, paste it from the clipboard. Note that after a Justify operation, this turns into unjustify. 

Ctrl+T: Check spelling. 

Ctrl+W: Find a word or phrase. At the prompt, use the cursor keys to go through previous search terms, or hit Ctrl+R to move into replace mode. Alternatively, you can hit Ctrl+T to go to a specific line.

Ctrl+C: Show current line number and file information. 

Ctrl+G: Get help; this provides information on navigating through files and common keyboard commands 
```

### Getting help in Linux

All Linux commands have manual pages. To access them, use “man” or “info” command. The manual page gives a detailed explanation of the command, all available options and sometimes, also provides examples. For example, to view the manual page for “ls” command:

```bash
man ls
```

Please explore manual pages of all the above commands for available options. 

### Linux text processing

*(a)* cut \
The cut command is a command line utility to cut a section from a file. Please see "man cut" for available options.

To cut a section of file use "-c" (characters)

```bash
cut -c1-10 SARS-CoV-2.fa
```

The option "-c1-10" will output first 10 characters from the input file. 

```
Few options: 
-c: cut based on character position 
-d: cut based on delimiter 
-f: field number 
```

We have a file named "human_viruses.txt" with some information including the names of the viruses, GenBank ids and genome length. These fields are separated by "|" symbol.

```bash
head human_viruses.txt
```

To get a list of the GenBank id,

```bash
cut -d "|" -f2 human_viruses.txt
```

*(b)* sort \
The sort command is used to sort the input content.

```
Few options: 
-t: field separator 
-n: numeric sort 
-k: sort with a key (field) 
-r: reverse sort 
-u: print unique entries 
```

```bash
sort -t "|" -nrk6 human_viruses.txt 
```

*(c)* grep \
grep searches the input for a given pattern.

```
Few options:
-A: after context
-B: before context
-C: before and after context
-c: count
-l: file with match
-i: ignore case
-o: only match
-v: invert match
-w: word match 
```

To get the list of all Hepatitis viruses from 'human_viruses.txt' file,

```bash
grep "Hepatitis" human_viruses.txt
```

*(d)* wc \
The command “wc” can be used in 2 ways, which counts lines, words or characters.

```bash
wc -l outbreak.csv
```

```bash
cat outbreak.csv | wc -l
```

*(e)* uniq \
The uniq command extracts unique lines from the input. It is usually used in combination with sort to count unique values in the input.

To get the list of countries that has had an outbreak in 2023:

```bash
cut -d, -f3 outbreak.csv | sort | uniq 
```

Other text processing commands worth looking at are: tr, rev, sed and paste.

### I/O control in Linux

When you run a command, the output is usually sent to standard output (stdout) ie. the terminal. However, we can redirect the standard output to a file using ">".

```bash
ls > list
cat list
```

The first command creates a new file called list with all the file names in the directory. If there exists a file already named "list", it is overwritten with the output of the command. Instead, we can append to a file using ">>" redirection. 

Another kind of output that is generated by programs is standard error. We must use “2>” to redirect it. 

ls foo 2> error

To redirect stdout and stderr to a file use "&>".

#### Pipes

Piping in Linux is a very powerful and efficient way to combine commands. Pipes (|) in Linux act as connecting links between commands. Pipe redirects output of the first command as an input to the next command. We can nest as many commands as we want using pipes. They ensure smooth running of the command flow and reduces the execution time. 

To print 10 smallest viruses,

```bash
sort -t"|" -nk6 human_viruses.txt | head -10
```

We will be working on other examples during the course, where we use pipes to combine more than two commands. 

### Process control

Some commands take time to finish the assigned job. For example, if you would like to compress a huge file with gzip command that takes a few minutes to finish running, you can run it in the background by appending the command with “&” (Another way is to suspend a command by pressing Ctrl+Z and typing “bg”). The completion of the task is indicated by “Done”.

```bash
gzip list &
```

We can get list of currently running jobs in the terminal by “jobs” command. This will give you all the background jobs running in the current terminal. If you want to see all the running processes in the system, use “top”. You can get user specific details in top using "-U" option. 

```bash
top
```

Few of the important columns in top output: 

* PID: Process Id, this is a unique number used to identify the process 
* COMMAND: Command Name 
* S: Process Status: The status of the task which can be one of: 
  + D = uninterruptible sleep 
  + R = running 
  + S = sleeping 
  + T = traced or stopped 
  + Z = zombie 

If you want to stop a running background job use “kill” command followed by the process id. 

kill 1234

This command kills the job with the process id 1234. As a user you can kill only your jobs. You do not have permission to run this command on the process ids of other users.

## Exercises

1. Open a new terminal and navigate into Exercises directory.
2. Extract first 15 lines from the file “NC_045512.2_gene_2.fa” and save the output into “output.fa”
3. How many fasta files are there in the directory?
4. Extract all header lines from the file "all.fa".
5. How many sequences are there in the file "all.fa"?
6. Copy the file "outbreak.csv" file into "Exercises" directory.
7. Get the list of countries that had at least two outbreaks in 2023.
8. Find the 99th line of the file "NC_045512.2_gene_1.fa" using only the 'tail' and 'head' command.
9. In the file "NC_078005.1_gene_1.fa", count the number of lines containing the sub-sequence "GGGG".
10. How do you stop a process with pid 5678?
11. Re-execute your previous command using a keyboard shortcut.
12. Create a new directory "Trial" and move "NC_012959.1_gene_4.fa" into the "Trial" directory.
13. Which is the command used to remove or delete file without confirmation prompt
14. __________ command is used to count the total number of lines, words and character in a file.
15. Which command would you use to extract 2nd, 5th, 7th column of a text file.
16. Extract first 10 lines of "outbreak.csv", sort them and save as "outbreak_1.csv". Extract first 20 lines of outbreak.csv, sort them and save as "outbreak_2.csv".
17. Extract common lines between the files "outbreak_1.csv" and "outbreak_2.csv" (use “comm” command, type "man comm" to get information) and save the output into a file named "output_1".
18. Which command would you use to find the word “pattern” from the file, “filename.txt”? Using that command, extract the “BioProject” information from the file SARS-CoV-2.gb.
19. Use the file SARS-CoV-2.gb to extract protein identifiers ("protein_id"). Remove the pattern "protein_id" from the output.
20. Which option with the command "rm" is required to remove a directory? 
