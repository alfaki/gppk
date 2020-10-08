/* gpplib.cpp (miscellaneous library routines) */

/******************************************************************************
 *  This code is part of GPPK (The Generalized Pooling Problem Kit).
 *
 *  Copyright (C) 2009, 2010 Mohammed Alfaki, Department of Informatics,
 *  University of Bergen, Bergen, Norway. All rights reserved. E-mail:
 *  <mohammeda@ii.uib.no>.
 *
 *  GPPK is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 ******************************************************************************/

#include "gppk.h"

/* Write formatted output to char string */
string get_str(const char *fmt, ...) {
	char *cstr;
	va_list arg;
	va_start(arg, fmt);
	int val = vasprintf(&cstr, fmt, arg);
	va_end(arg);
	string str(cstr);
	free(cstr);
	return str;
}

/* print error message and terminate processing */
void xerror(const char *fmt, ...) {
	va_list arg;
	va_start(arg, fmt);
	vprintf(fmt, arg);
	va_end(arg);
	printf("\n");
	longjmp(buf, 1);
	/* no return */
}

/* return file name for a full path of file */
string file_name(const char *file) {
	string str = file, fname, fpath;
	size_t size = str.find_last_of("/\\");
	fpath = str.substr(0,size);
	fname = str.substr(size+1);
	size = str.substr(size+1).find_last_of(".");
	fname = fname.substr(0, size);
	return fname;
}

/* return the current date and time char string */
char *cdate() {
	time_t rawtime;
	struct tm *timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	return asctime(timeinfo);
}

/* convert time in seconds to hh:mm:ss */
string clocktime(double ttime) {
	char *ctime;
	ctime = new char[STRING_LENGTH];
	int time, hr, mn, sc, ms;
	time = (int)ttime;
	hr = time/3600; time = time%3600;
	mn = time/60; sc = time%60;
	ms = (int)(100*(ttime-(double)time));
	sprintf(ctime, "%02d:%02d:%02d:%02d", hr, mn, sc, ms);
	string str(ctime);
	delete[] ctime;
	return str;
}

/* compute n choose k */
int n_choose_k(int n, int k) {
	if (k <= 0 || k >= n)
		return 1;
	else
		return n_choose_k(n-1, k-1)+n_choose_k(n-1, k);
}

/******************************************************************************
 *  NAME
 *  str2int - convert character string to value of int type
 *
 *  SYNOPSIS
 *  #include "gppk.h"
 *  int str2int(const char *str, int *val);
 *
 *  DESCRIPTION
 *  The routine str2int converts the character string str to a value of
 *  integer type and stores the value into location, which the parameter
 *  val points to (in the case of error content of this location is not
 *  changed).
 *
 *  RETURNS
 *  The routine returns one of the following error codes:
 *  0 - no error;
 *  1 - value out of range;
 *  2 - character string is syntactically incorrect.
 ******************************************************************************/

int str2int(const char *str, int *_val) {
	int d, k, s, val = 0;
	/* scan optional sign */
	if (str[0] == '+')
		s = +1, k = 1;
	else if (str[0] == '-')
		s = -1, k = 1;
	else
		s = +1, k = 0;
	/* check for the first digit */
	if (!isdigit((unsigned char)str[k]))
		return 2;
	/* scan digits */
	while (isdigit((unsigned char)str[k])) {
		d = str[k++] - '0';
		if (s > 0) {
			if (val > INT_MAX / 10)
				return 1;
			val *= 10;
			if (val > INT_MAX - d)
				return 1;
			val += d;
		}
		else {
			if (val < INT_MIN / 10) return 1;
			val *= 10;
			if (val < INT_MIN + d) return 1;
			val -= d;
		}
	}
	/* check for terminator */
	if (str[k] != '\0')
		return 2;
	/* conversion has been done */
	*_val = val;
	return 0;
}

/******************************************************************************
 *  NAME
 *  str2num - convert character string to value of double type
 *
 *  SYNOPSIS
 *  #include ".h"
 *  int str2num(const char *str, double *val);
 *
 *  DESCRIPTION
 *  The routine str2num converts the character string str to a value of
 *  double type and stores the value into location, which the parameter
 *  val points to (in the case of error content of this location is not
 *  changed).
 *
 *  RETURNS
 *  The routine returns one of the following error codes:
 *  0 - no error;
 *  1 - value out of range;
 *  2 - character string is syntactically incorrect.
 ******************************************************************************/

int str2num(const char *str, double *_val) {
	int k;
	double val;
	/* scan optional sign */
	k = (str[0] == '+' || str[0] == '-' ? 1 : 0);
	/* check for decimal point */
	if (str[k] == '.') {
		k++;
		/* a digit should follow it */
		if (!isdigit((unsigned char)str[k]))
			return 2;
		k++;
		goto frac;
	}
	/* integer part should start with a digit */
	if (!isdigit((unsigned char)str[k]))
		return 2;
	/* scan integer part */
	while (isdigit((unsigned char)str[k]))
		k++;
	/* check for decimal point */
	if (str[k] == '.')
		k++;
	frac: /* scan optional fraction part */
	while (isdigit((unsigned char)str[k]))
		k++;
	/* check for decimal exponent */
	if (str[k] == 'E' || str[k] == 'e') {
		k++;
		/* scan optional sign */
		if (str[k] == '+' || str[k] == '-')
			k++;
		/* a digit should follow E, E+ or E- */
		if (!isdigit((unsigned char)str[k]))
			return 2;
	}
	/* scan optional exponent part */
	while (isdigit((unsigned char)str[k]))
		k++;
	/* check for terminator */
	if (str[k] != '\0')
		return 2;
	/* perform conversion */
	{
		char *endptr;
		val = strtod(str, &endptr);
		if (*endptr != '\0')
			return 2;
	}
	/* check for overflow */
	if (!(-DBL_MAX <= val && val <= +DBL_MAX))
		return 1;
	/* check for underflow */
	if (-DBL_MIN < val && val < +DBL_MIN)
		val = 0.0;
	/* conversion has been done */
	*_val = val;
	return 0;
}

/* eof */
