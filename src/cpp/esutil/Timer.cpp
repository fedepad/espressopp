#include "Timer.hpp"

#ifdef HAVE_SYS_RESOURCE_H

#include <sys/resource.h>

real esutil::UserTimer::getCurrentTime() const {
  struct rusage rus;
  getrusage(RUSAGE_SELF, &rus);
  return rus.ru_utime.tv_sec + 1.e-6*rus.ru_utime.tv_usec;
}
#else

// we do not have getrusage
real esutil::UserTimer::getCurrentTime() const { return 0; }
#endif

#ifdef HAVE_BOOST_MPI

real esutil::WallTimer::getCurrentTime() const { return timer.elapsed(); }

#elif defined(HAVE_SYS_TIME_H)

#include <sys/time.h>

real esutil::WallTimer::getCurrentTime() const {
  struct timeval tp;
  struct timezone tzp;
  gettimeofday (&tp, &tzp);
  return tp.tv_sec + 1.e-6*tp.tv_usec;
}

#elif defined(HAVE_TIME_H)

#include <time.h>

real esutil::WallTimer::getCurrentTime() const { return time(0); }

#else

// we do not have gettimeofday
real esutil::WallTimer::getCurrentTime() const { return 0; }

#endif
