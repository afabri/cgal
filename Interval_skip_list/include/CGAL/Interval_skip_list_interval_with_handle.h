// Copyright (c) 2003 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_INTERVAL_SKIP_LIST_INTERVAL_WITH_HANDLE_H
#define CGAL_INTERVAL_SKIP_LIST_INTERVAL_WITH_HANDLE_H

#include <CGAL/license/Interval_skip_list.h>
#include <CGAL/Interval_skip_list_interval.h>

#include <CGAL/basic.h>
#include <cstdlib>
#include <iostream>


namespace CGAL {

  template <class Value_, class H>
  class Interval_skip_list_interval_with_handle : public Interval_skip_list_interval<Value_>
  {
  public:
    typedef Value_ Value;
    typedef Interval_skip_list_interval<Value_> Base;

  private:
    H handle_;
    bool lbound_;
    bool rbound_;
    Value inf_;
    Value sup_;
  public:

    Interval_skip_list_interval_with_handle()
      : Base()
    {}

    Interval_skip_list_interval_with_handle(H handle,
                                            const Value& inf_,
                                            const Value& sup_,
                                            bool lb = true,
                                            bool rb = true)
      : Base(inf_, sup_,lb,rb), handle_(handle)
    {}

    H handle() const { return handle_;}

    bool operator==(const Interval_skip_list_interval_with_handle& I) const
    {
      return ( (handle() == I.handle()) &&
               (this->inf() == I.inf()) && (this->sup() == I.sup()) &&
	       (this->inf_closed() == I.inf_closed()) && (this->sup_closed() == I.sup_closed()) );
    }

    bool operator!=(const Interval_skip_list_interval_with_handle& I) const
    {
      return ! (*this == I);
    }
  };
} // namespace CGAL

#endif // CGAL_INTERVAL_SKIP_LIST_INTERVAL_WITH_HANDLE_H
