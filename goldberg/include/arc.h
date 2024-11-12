// Hey emacs! This is -*- c++ -*-
// $Id: arc.h,v 1.5 1997/07/01 15:25:44 mslevine Exp $

#ifndef ARC_H
#define ARC_H

/***********************************************************************/
/*                                                                     */
/*         CLASS arc_basic_data                                        */
/*                                                                     */
/***********************************************************************/
class graph;
class arc_basic_data
{
  node *hd;

public:
  void setHead(node *v) {hd=v;}
  node *head()   {return hd;}
};

#endif


