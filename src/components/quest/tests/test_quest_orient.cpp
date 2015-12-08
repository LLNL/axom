/*
 * $Id$
 */

/*!
 *******************************************************************************
 * \file
 *
 * \date Dec 18, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

#include "gtest/gtest.h"

#include "quest/Orientation.hpp"
#include "quest/Point.hpp"
#include "quest/Triangle.hpp"


TEST( quest_orientation, orient3D )
{
   // STEP 0: Setup triangle ABC in 3D
   quest::Point< double,3 > A =
           quest::Point< double,3 >::make_point( 0.0, 0.0, 0.0 );
   quest::Point< double,3 > B =
           quest::Point< double,3 >::make_point( 1.5, 1.5, 0.0 );
   quest::Point< double,3 > C =
           quest::Point< double,3 >::make_point( 2.5, 0.0, 0.0 );
   quest::Triangle< double,3 > tri(A,B,C);

   // STEP 1: Setup test points, q0, q1, q2 => boundary, positive, negative
   quest::Point< double,3 > q0=
           quest::Point< double,3 >::make_point(1.5,0.5,0.0);
   quest::Point< double,3 > q1=
           quest::Point< double,3 >::make_point(1.5,0.5,0.5);
   quest::Point< double,3 > q2 =
           quest::Point< double,3 >::make_point(1.5,0.5,-0.5);

   // STEP 2: test orientation
   int orient = quest::orientation( q0, tri );
   EXPECT_EQ( quest::ON_BOUNDARY,orient);

   orient = quest::orientation( q1, tri );
   EXPECT_EQ( quest::ON_NEGATIVE_SIDE, orient );

   orient = quest::orientation( q2, tri );
   EXPECT_EQ( quest::ON_POSITIVE_SIDE, orient );
}


