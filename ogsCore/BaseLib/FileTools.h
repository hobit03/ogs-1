/**
 * \file FileTools.h
 * 26/4/2010 LB Initial implementation
 *
 */


#ifndef FILETOOLS_H
#define FILETOOLS_H

// ** INCLUDES **
#include <sys/stat.h>

namespace BaseLib {
/**
 * Returns true if given file exists. From http://www.techbytes.ca/techbyte103.html
 */
static bool IsFileExisting(const std::string &strFilename)
{
	struct stat stFileInfo;
	bool blnReturn;
	int intStat;

	// Attempt to get the file attributes
	intStat = stat(strFilename.c_str(),&stFileInfo);

	if(intStat == 0)
	{
		// We were able to get the file attributes
		// so the file obviously exists.
		blnReturn = true;
	}
	else
	{
		// We were not able to get the file attributes.
		// This may mean that we don't have permission to
		// access the folder which contains this file. If you
		// need to do that level of checking, lookup the
		// return values of stat which will give you
		// more details on why stat failed.
		blnReturn = false;
	}

	return(blnReturn);
}


static std::string getFileDirecotryPath(const std::string &file_path)
{
	size_t indexChWin, indexChLinux;
	indexChWin = indexChLinux = 0;
	indexChWin = file_path.find_last_of('\\');
	indexChLinux = file_path.find_last_of('/');
	//
	std::string dir_path;
	if(indexChWin != std::string::npos)
		dir_path = file_path.substr(0,indexChWin) + "\\";
	else if(indexChLinux != std::string::npos)
		dir_path = file_path.substr(0,indexChLinux) + "/";

	return dir_path;
}

static std::string getFileBaseName(const std::string &file_path)
{
	size_t indexChWin, indexChLinux;
	indexChWin = indexChLinux = 0;
	indexChWin = file_path.find_last_of('\\');
	indexChLinux = file_path.find_last_of('/');
	//
	std::string dir_path;
	if(indexChWin != std::string::npos)
		dir_path = file_path.substr(indexChWin, file_path.length());
	else if(indexChLinux != std::string::npos)
		dir_path = file_path.substr(indexChLinux, file_path.length());

	return dir_path;
}

} // end namespace BaseLib

#endif // FILETOOLS_H
