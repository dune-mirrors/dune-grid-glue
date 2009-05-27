// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    timer.hh
 *  Version:     1.0
 *  Created on:  Jan 24, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: Class for measuring runtimes.
 *  subversion:  $Id$
 *
 */
/**
 * @file timer.hh
 * @brief simple support for measuring usertime
 */

#ifndef TIMER_HH_
#define TIMER_HH_

// headers for getrusage (2)
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

// headers for stderror (3)
#include <string.h>
// access to errno
#include <errno.h>
#include <iostream>


/**
 * @class UserTimeMeasure
 * @brief Class for measuring runtimes.
 *
 * Only measures time that the cpu is in use.
 */
template<typename OS>
class UserTimeMeasure
{
public:

  UserTimeMeasure(OS &os, const char* name) : _os(os), _name(name)
  {}

  void setName(const char* name)
  {
    this->_name = name;
  }

  //! reset timer
  void reset()
  {
    struct rusage ru;
    if (getrusage(RUSAGE_SELF, &ru))
      this->_os << "Error in '" << this->_name << "': " << strerror(errno) << std::endl;
    this->_timer = ru.ru_utime;
  }

  //! get elapsed user-time in seconds
  double getTime()
  {
    struct rusage ru;
    if (getrusage(RUSAGE_SELF, &ru))
      this->_os << "Error in '" << this->_name << "': " << strerror(errno) << std::endl;
    //		return 1.0 * (ru.ru_utime.tv_sec - this->_timer.tv_sec) + 1E-3*(ru.ru_utime.tv_usec - this->_timer.tv_usec);
    return 1.0 * (ru.ru_utime.tv_sec - this->_timer.tv_sec) + 1E-6*(ru.ru_utime.tv_usec - this->_timer.tv_usec);
  }

  //! print elapsed user-time in seconds
  void print()
  {
    this->_os << "Timer '" << this->_name << "' elapsed time: " << this->getTime() << std::endl;
  }

private:

  struct timeval _timer;
  OS& _os;
  std::string _name;
};

#endif // TIMER_HH_
