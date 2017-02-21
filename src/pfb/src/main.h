
#ifndef __main_H__
#define __main_H__

#include <stdio.h>
#include <stdlib.h>

#include <string.h>		/* for memset(), strncopy(), memcpy(), strerror() */
#include <sys/types.h>	/* for open() */
#include <sys/stat.h>   /* for open() */
#include <fcntl.h>		/* for open() */
#include <unistd.h>		/* for close() */

#include <float.h>		/* for FLT_MAX */
#include <getopt.h>		/* for option parsing */
#include <assert.h>		/* for assert() */
#include <errno.h>		/* for errno*/
#include <signal.h>		/* for signal-hadnling */

#include <sys/time.h>	/* for gettimeofday() */