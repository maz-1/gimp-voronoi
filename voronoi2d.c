/*
 * @(#) $Id: voronoi2d.c,v 1.15 2004/05/08 19:13:22 yeti Exp $
 * Generate 2D Voronoi diagrams in linear time.
 *
 * Copyright (C) 2001-2002,2004 Yeti (David Necas) <yeti@physics.muni.cz>.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 */
#define _GNU_SOURCE
#include <stdlib.h>
#include <math.h>
#include <glib.h>
#include "voronoi.h"

/* TODO: */
gdouble radial_factor;

static inline gdouble
angle(const VoronoiCoords *const r)
{
  return atan2(r->y, r->x);
}

static gint
vobj_angle_compare(gconstpointer x, gconstpointer y)
{
  gdouble xangle, yangle;

  xangle = ((VoronoiObject*)x)->angle;
  yangle = ((VoronoiObject*)y)->angle;

  if (xangle < yangle)
    return -1;
  if (xangle > yangle)
    return 1;
  return 0;
}

/* XXX: owner->ne requirements: NONE */
gdouble
value_flat(const VoronoiCoords *const point,
           const VoronoiObject *const owner)
{
  gdouble r;

  r = owner->random;

  return r;
}

/* XXX: owner->ne requirements: NONE */
gdouble
value_radial(const VoronoiCoords *const point,
             const VoronoiObject *const owner)
{
  VoronoiCoords dist;
  gdouble r;

  dist = coords_minus(point, &owner->pos);
  r = radial_factor*sqrt(DOTPROD_SS(dist, dist));

  return r;
}

/* XXX: owner->ne requirements: cyclic, neighbourized, segment angles */
gdouble
value_segmented(const VoronoiCoords *const point,
                const VoronoiObject *const owner)
{
  VoronoiCoords dist;
  gdouble r, phi;
  GSList *ne;

  ne = owner->ne;
  dist = coords_minus(point, &owner->pos);
  phi = angle(&dist);

  while ((phi >= VOBJ(ne)->angle)
         + (phi < VOBJ(ne->next)->angle)
         + (VOBJ(ne)->angle > VOBJ(ne->next)->angle) < 2)
    ne = ne->next;

  r = 2*DOTPROD_SS(dist, VOBJ(ne)->rel.v)/VOBJ(ne)->rel.d;

  return r;
}

/* XXX: owner->ne requirements: neighbourized */
gdouble
value_border(const VoronoiCoords *const point,
             const VoronoiObject *const owner)
{
  VoronoiCoords dist;
  gdouble r, r_min;
  GSList *ne;

  dist = coords_minus(point, &owner->pos);
  r_min = HUGE_VAL;

  for (ne = owner->ne; ne != NULL; ne = ne->next) {
    r = fabs(VOBJ(ne)->rel.d/2 - DOTPROD_SS(dist, VOBJ(ne)->rel.v))
        /sqrt(VOBJ(ne)->rel.d);
    if (r < r_min)
      r_min = r;
    if (ne->next == owner->ne)
      break;
  }

  r = 1 - 2*r_min*radial_factor;

  return r;
}

/* XXX: owner->ne requirements: NONE */
gdouble
value_second(const VoronoiCoords *const point,
             const VoronoiObject *const owner)
{
  VoronoiCoords dist;
  gdouble r, r_min;
  GSList *ne;

  r_min = HUGE_VAL;

  for (ne = owner->ne; ne != NULL; ne = ne->next) {
    dist = coords_minus(point, &VOBJ(ne)->pos);
    r = DOTPROD_SS(dist, dist);
    if (r < r_min)
      r_min = r;
    if (ne->next == owner->ne)
      break;
  }

  r = 1 - sqrt(r_min)*radial_factor;

  return r;
}

/* XXX: owner->ne requirements: NONE */
gdouble
value_sprod(const VoronoiCoords *const point,
            const VoronoiObject *const owner)
{
  VoronoiCoords dist, dist_min;
  gdouble r, r_min;
  GSList *ne;

  r_min = HUGE_VAL;
  dist_min.x = dist_min.y = 0.0;

  for (ne = owner->ne; ne != NULL; ne = ne->next) {
    dist = coords_minus(point, &VOBJ(ne)->pos);
    r = DOTPROD_SS(dist, dist);
    if (r < r_min) {
      r_min = r;
      dist_min = dist;
    }
    if (ne->next == owner->ne)
      break;
  }

  dist = coords_minus(&owner->pos, point);
  r = (DOTPROD_SS(dist, dist_min)*radial_factor*radial_factor + 0.5)/1.5;

  return r;
}

/* compute segment angles
   more precisely, VOBJ(ne)->angle will be set to start angle for segment
   from ne to ne->next (so end angle is in ne->next)
   XXX: ne0 requirements: cyclic and neighbourized */
void
compute_segment_angles(GSList *ne0)
{
  GSList *ne;
  VoronoiObject *p, *q;
  VoronoiCoords z;

  ne = ne0;
  do {
    p = VOBJ(ne);
    q = VOBJ(ne->next);
    z.x = p->rel.d * q->rel.v.y - q->rel.d * p->rel.v.y;
    z.y = q->rel.d * p->rel.v.x - p->rel.d * q->rel.v.x;
    q->angle = angle(&z);
    ne = g_slist_next(ne);
  } while (ne != ne0);
}

/* compute angles from rel.v relative coordinates
   XXX: ne0 requirements: neighbourized */
void
compute_straight_angles(GSList *ne0)
{
  GSList *ne;
  VoronoiObject *p;

  for (ne = ne0; ne; ne = g_slist_next(ne)) {
    p = VOBJ(ne);
    p->angle = angle(&p->rel.v);
    if (ne->next == ne0)
      return;
  }
}

/* compute relative positions and norms to center center
   XXX: ne0 requirements: NONE */
void
neighbourize(GSList *ne0, const VoronoiCoords *const center)
{
  GSList *ne;

  for (ne = ne0; ne; ne = g_slist_next(ne)) {
    VoronoiObject *p = VOBJ(ne);

    p->rel.v = coords_minus(&p->pos, center);
    p->rel.d = DOTPROD_SS(p->rel.v, p->rel.v);
    if (ne->next == ne0)
      return;
  }
}

/* return intersection time t for intersection of lines:

           r = linevec*t + start
     |r - a| = |r - b|
 */
static inline gdouble
intersection_time(const VoronoiCoords *const a, const VoronoiCoords *const b,
                  const VoronoiCoords *const linevec,
                  const VoronoiCoords *const start)
{
  VoronoiCoords p, q;
  gdouble s;

  /* line dividing a-neighbourhood and b-neighbourhood */
  q = coords_minus(b, a);
  p = coords_plus(b, a);

  /* FIXME: this is numerically unstable.  I didn't observed any problems even
     for sigma = 0.00001, nevertheless it is. (the probability something odd
     will happen fortunately decreases fast enough with increasing oddity of
     the thing)  almost zero denominator fortunately doesn't matter much, as
     long as numerator is not almost zero too.  then we are trying to find
     intersection point of two almost parallel lines, which is moreover quite
     near, so we want to know its position precisely.  I don't know how to do
     it with the least loss of precession.  sigh. */
  s = DOTPROD_SP(q, linevec);
  if (fabs(s) < 1e-14)
    s = 1e-14; /* better than nothing */
  return (DOTPROD_SS(q, p)/2 - DOTPROD_SP(q, start))/s;
}

/* being in point start owned by owner (XXX: this condition MUST be true)
   we want to get to point end and know our new owner
   returns the new owner; in addition, when next_safe is not NULL it stores
   there number of times we can repeat move along (end - start) vector still
   remaining in the new owner */
VoronoiObject*
move_along_line(const VoronoiObject *const owner,
                const VoronoiCoords *const start,
                const VoronoiCoords *const end, gint *next_safe)
{
  VoronoiCoords linevec;
  VoronoiObject *ow;
  GSList *ne, *nearest = NULL;
  gdouble t, t_min, t_back;

  ow = (VoronoiObject*)owner;
  linevec = coords_minus(end, start);
  t_back = 0;
  /* XXX: start must be owned by owner, or else strange things will happen */
  for ( ; ; ) {
    t_min = HUGE_VAL;
    ne = ow->ne;
    do {
      /* find intersection with border line between ow and ne
         FIXME: there apparently exist values t > t_back && t_back > t */
      t = intersection_time(&ow->pos, &VOBJ(ne)->pos, &linevec, start);
      if (t - t_back >= EPS && t < t_min) {
        t_min = t;
        nearest = ne;
      }
      ne = ne->next;
    } while (ne != ow->ne);

    /* no intersection inside the abscissa? then we are finished and can
       compute how many steps the same direction will remain in ow's
       neighbourhood */
    if (t_min > 1) {
      if (next_safe == NULL)
        return ow;
      if (t_min == HUGE_VAL)
        *next_safe = G_MAXINT;
      else
        *next_safe = floor(t_min) - 1;
      return ow;
    }

    /* otherwise nearest intersection determines a new owner */
    ow = VOBJ(nearest);
    t_back = t_min; /* time value showing we are going back */
  }
}

/* find and return owner of point point
   NB: this is crude and should not be used for anything else than initial
   grip, use move_along_line() then
   N: works for both cyclic and noncyclic ne-> */
VoronoiObject*
find_owner(GSList **squares,
           const gint wsq,
           const gint hsq,
           const VoronoiCoords *const point)
{
  GSList *ne;
  VoronoiObject *owner = NULL;
  VoronoiCoords dist;
  gint jx, jy;
  gint ix, iy;
  gdouble norm_min;

  jx = floor(point->x) + SQBORDER;
  jy = floor(point->y) + SQBORDER;

  /* scan the 25-neighbourhood
     XXX: I thought 9-neighbourhood was enough, but apparently it wasn't */
  norm_min = HUGE_VAL;
  for (ix = -2; ix <= 2; ix++) {
    if (jx + ix < 0 || jx + ix >= wsq + 2*SQBORDER)
      continue;
    for (iy = -2; iy <= 2; iy++) {
      if (jy + iy < 0 || jy + iy >= hsq + 2*SQBORDER)
        continue;
      for (ne = squares[(iy+jy)*(wsq + 2*SQBORDER) + ix+jx];
           ne != NULL;
           ne = ne->next) {
        dist = coords_minus(&VOBJ(ne)->pos, point);
        if (DOTPROD_SS(dist, dist) < norm_min) {
          norm_min = DOTPROD_SS(dist, dist);
          owner = VOBJ(ne);
        }
        if (ne->next == squares[(iy+jy)*(wsq + 2*SQBORDER) + ix+jx])
          break;
      }
    }
  }

  return owner;
}

/* return true iff point z (given as VoronoiLine) is shadowed by points a and b
   (XXX: all coordiantes are relative) */
static inline gboolean
in_shadow(const VoronoiLine *const a, const VoronoiLine *const b,
          const VoronoiCoords *const z)
{
  VoronoiCoords r, oa, ob, rz;
  gdouble s;

  /* artifical fix for periodic grids, because in Real World This Just Does Not
   * Happen
   * also mitigates the s == 0 case below, as the offending point would
   * be probably removed here */
  if (DOTPROD_SP(a->v, z) > 1.01*a->d
      && fabs(a->v.x * z->y - z->x * a->v.y) < 1e-12)
    return TRUE;
  if (DOTPROD_SP(b->v, z) > 1.01*b->d
      && fabs(b->v.x * z->y - z->x * b->v.y) < 1e-12)
    return TRUE;

  s = 2*(a->v.x * b->v.y - b->v.x * a->v.y);
  /* FIXME: what to do when s == 0 (or very near)??? */
  r.x = (a->d * b->v.y - b->d * a->v.y)/s;
  r.y = (b->d * a->v.x - a->d * b->v.x)/s;
  oa.x = -a->v.y;
  oa.y = a->v.x;
  ob.x = -b->v.y;
  ob.y = b->v.x;
  rz = coords_minus(z, &r);
  return (DOTPROD_SS(rz, rz) > DOTPROD_SS(r, r)
          && DOTPROD_PS(z, oa)*DOTPROD_SS(b->v, oa) > 0
          && DOTPROD_PS(z, ob)*DOTPROD_SS(a->v, ob) > 0);
}

static GSList*
extract_neighbourhood(GSList **squares,
                      gint wsq,
                      gint hsq,
                      VoronoiObject *p)
{
  GSList *ne = NULL;
  gint jx, jy;
  gint ix, iy;
  gint xwsq, xhsq;

  xwsq = wsq + 2*SQBORDER;
  xhsq = hsq + 2*SQBORDER;

  jx = floor(p->pos.x) + SQBORDER;
  jy = floor(p->pos.y) + SQBORDER;

  /* construct the 37-neighbourhood list */
  for (ix = -3; ix <= 3; ix++) {
    if (jx + ix < 0 || jx + ix >= xwsq)
      continue;
    for (iy = -3; iy <= 3; iy++) {
      if (abs(ix) == 3 && abs(iy) == 3)
        continue;
      if (jy + iy < 0 || jy + iy >= xhsq)
        continue;
      ne = g_slist_concat(g_slist_copy(squares[(iy+jy)*xwsq + ix+jx]), ne);
      if (ix == 0 && iy == 0)
        ne = g_slist_remove(ne, p);
    }
  }

  g_assert(ne != NULL);

  /* compute relative coordinates and angles */
  neighbourize(ne, &p->pos);
  compute_straight_angles(ne);

  return ne;
}

static GSList*
shadow_filter(GSList *ne)
{
  GSList *ne1, *ne2;
  gint notremoved;
  gint len;

  if (ne == NULL)
    return ne;

  /* make the list cyclic if it isn't already
   * (we have to unlink elements ourself then) */
  len = 1;
  for (ne2 = ne; ne2->next && ne2->next != ne; ne2 = g_slist_next(ne2))
    len++;
  if (len < 3)
    return ne;
  ne2->next = ne;

  /* remove objects shadowed by their ancestors and successors
     XXX: in non-degenerate case this is O(n*log(n)), but can be O(n*n) */
  ne1 = ne;
  notremoved = 0;
  do {
    ne2 = ne1->next;
    if (in_shadow(&VOBJ(ne1)->rel,
                  &VOBJ(ne2->next)->rel,
                  &VOBJ(ne2)->rel.v)) {
      ne1->next = ne2->next;
      g_slist_free_1(ne2);
      notremoved = 0;
      len--;
    }
    else {
      ne1 = ne2;
      notremoved++;
    }
  } while (notremoved < len && len > 2);

  return ne1; /* return cyclic list */
}

/* object-removal code, not used as the square side seems to be large enough
 * to make the number of objects artifically inserted into otherwise empty
 * squares really negligible => there's no need to simulate empty squares
 * by removing the objects later */
#if 0
/* duplicate possibly cyclic list as either cyclic or noncyclic */
static inline GSList*
duplicate_list(GSList *list,
               gboolean create_cyclic)
{
  GSList *l, *copy;

  if (!list)
    return NULL;

  copy = NULL;
  for (l = list; l->next && l->next != l; l = g_slist_next(l))
    copy = g_slist_prepend(copy, l->data);
  if (l != list)
    copy = g_slist_prepend(copy, l->data);

  l = copy;
  copy = g_slist_reverse(copy);
  if (create_cyclic)
    l->next = copy;

  return copy;
}

/* compute the number of squares with just one object
 * each square MUST contain some */
static gint
count_1squares(GSList **squares,
               const gint wsq,
               const gint hsq)
{
  gint i, count;
  gint xwsq, xhsq;

  xwsq = wsq + 2*SQBORDER;
  xhsq = hsq + 2*SQBORDER;

  count = 0;
  for (i = 0; i < xwsq*xhsq; i++) {
    if (!squares[i]->next)
      count++;
  }

  return count;
}

static void
remove_object(VoronoiObject *obj)
{
  GSList *necopy, *ne, *nene, *ne2;

  g_return_if_fail(obj->ne);
  /* neighbours of the removed object */
  necopy = duplicate_list(obj->ne, FALSE);
  for (ne = necopy; ne; ne = g_slist_next(ne)) {
    VoronoiObject *neighbour = VOBJ(ne);

    /* join both neighbourhoods */
    nene = g_slist_concat(duplicate_list(neighbour->ne, FALSE),
                          g_slist_copy(necopy));
    /* remove object and neighbour itself */
    nene = g_slist_remove_all(nene, obj);
    nene = g_slist_remove_all(nene, neighbour);
    /* neighbourize */
    neighbourize(nene, &neighbour->pos);
    compute_straight_angles(nene);
    nene = g_slist_sort(nene, &vobj_angle_compare);
    /* remove duplicities */
    for (ne2 = nene; ne2 && ne2->next; ne2 = g_slist_next(ne2)) {
      if (ne2->data == ne2->next->data) {
        GSList *nex;

        nex = ne2->next;
        ne2->next = nex->next;
        g_slist_free_1(nex);
      }
    }
    /* shadow filter and replace the neighbourhood */
    nene = shadow_filter(nene);
    g_slist_free(neighbour->ne);
    neighbour->ne = nene;
  }
  g_slist_free(necopy);
}

void
remove_objects(GSList **squares,
               const gint wsq,
               const gint hsq,
               gboolean htil,
               gboolean vtil)
{
  gdouble n1;
  gint nremove;
  gint xwsq, xhsq;

  xwsq = wsq + 2*SQBORDER;
  xhsq = hsq + 2*SQBORDER;

  n1 = count_1squares(squares, wsq, hsq);
  nremove = n1*n1/(xwsq*xhsq);
  fprintf(stderr, "n1 = %.0f, N = %d, nremove = %d\n",
          n1, xwsq*xhsq, nremove);
}
#endif

void
find_voronoi_neighbours_iter(GSList **squares,
                             const gint wsq,
                             const gint hsq,
                             const gint iter)
{
  GSList *this;

  for (this = squares[iter]; this; this = g_slist_next(this)) {
    VoronoiObject *obj = VOBJ(this);

    obj->ne = extract_neighbourhood(squares, wsq, hsq, obj);
    obj->ne = g_slist_sort(obj->ne, &vobj_angle_compare);
    obj->ne = shadow_filter(obj->ne);
  }
}

void
find_voronoi_neighbours(GSList **squares,
                        const gint wsq,
                        const gint hsq)
{
  gint xwsq, xhsq, i;

  xwsq = wsq + 2*SQBORDER;
  xhsq = hsq + 2*SQBORDER;

  for (i = 0; i < xwsq*xhsq; i++)
    find_voronoi_neighbours_iter(squares, wsq, hsq, i);
}

void
destroy_neighbourhoods(GSList **squares,
                       const gint wsq,
                       const gint hsq)
{
  gint xwsq, xhsq, i;
  GSList *this, *ne;

  xwsq = wsq + 2*SQBORDER;
  xhsq = hsq + 2*SQBORDER;

  for (i = 0; i < xwsq*xhsq; i++) {
    for (this = squares[i]; this; this = g_slist_next(this)) {
      ne = VOBJ(this)->ne->next;
      VOBJ(this)->ne->next = NULL; /* break cycle */
      g_slist_free(ne);
    }
  }
}

