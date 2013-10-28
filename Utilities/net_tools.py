'''
network.py

Synopsis:
	Contains useful functions relating to network connectivity
	
Requires:
	pexpect module
	http://www.noah.org/wiki/pexpect
	
Comments:
	This is super useful! Not so secure though...
'''


import pexpect
import sys



def secure_copy (FILE, REMOTE_LOC, REMOTE_HOST, USER, PASSWORD, FOLDER=False):

	'''
	Copies a file from your local location to a remote location
	
	Arguments:
		FILE							path/to/file to copy
		REMOTE_LOC					Remote path
		REMOTE_HOST					Remote host
		USER							Remote username
		PASSWORD						Password (note irony in 'secure')
		FOLDER = True, False			Are you copying a folder? If so will use -r flag, , default False
	
	
	Returns:
		N/A
		
	Example:
		secure_copy ( 'file.txt', '/home/jm8g08/', 'cosmos.phys.soton.ac.uk', 'jmatthews', 'mypasswd')
	'''
	
	if Folder: 
		scp = "scp -r -oPubKeyAuthentication=no"
	else:
		scp= "scp -oPubKeyAuthentication=no"

	COMMAND = "%s %s %s@%s:%s" % (scp, FILE, USER, HOST, REMOTE_FILE)

	child = pexpect.spawn ( COMMAND )
	child.expect ('password:')
	child.sendline ( PASS )
	child.expect ( pexpect.EOF )



def secure_get (REMOTE_FILE, REMOTE_HOST, USER, PASSWORD, FOLDER=False, LOC="." ):

	'''
	Gets a file from a remote location to your local location 
	
	Arguments:
		REMOTE_FILE					path/to/file to copy on remote machine
		REMOTE_HOST					Remote host
		USER							Remote username
		PASSWORD						Password (note irony in 'secure')
		LOC							Local path, default current dir
		FOLDER = True, False			Are you copying a folder? If so will use -r flag, default False
	
	
	Returns:
		N/A
		
	Example:
		secure_get ( '/home/jm8g08/file.txt', 'cosmos.phys.soton.ac.uk', 'jmatthews', 'mypasswd')
	'''
	
	if Folder: 
		scp = "scp -r -oPubKeyAuthentication=no"
	else:
		scp= "scp -oPubKeyAuthentication=no"

	COMMAND = "%s %s@%s:%s %s" % (scp, USER, HOST, REMOTE_FILE, LOC)

	child = pexpect.spawn ( COMMAND )
	child.expect ('password:')
	child.sendline ( PASS )
	child.expect ( pexpect.EOF )
	
	
