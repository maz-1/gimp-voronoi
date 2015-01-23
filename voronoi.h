/* @(#) $Id: voronoi.h,v 1.15 2004/05/08 19:13:22 yeti Exp $ */
#ifndef VORONOI_H
#define VORONOI_H

#include <glib.h>

#include <stdio.h>

/* just some small but definitely nonzero number (it should be several orders
   larger than machine epsilon, at least) */
#define EPS 0.0000001
/* how larger the squarized grid should be (measured in squares) */
#define SQBORDER 2

extern GRand *rng;

/* grid type
   (more precisely, neighbourhood type---square grid is really square, but to
   get triangular neighbourhood the grid itself must be hexagonal and vice
   versa; random grid is really random) */
typedef enum {
  VORONOI_GRID_RANDOM = -1,
  VORONOI_GRID_SQUARE = 0,
  VORONOI_GRID_RHOMB = 1,
  VORONOI_GRID_TRIANGLE_H = 2,
  VORONOI_GRID_TRIANGLE_V = 3,
  VORONOI_GRID_HEXAGON_H = 4,
  VORONOI_GRID_HEXAGON_V = 5,
} VoronoiGridType;

typedef enum {
  VORONOI_VALUE_FLAT = 0,
  VORONOI_VALUE_RADIAL = 1,
  VORONOI_VALUE_SEGMENTED = 2,
  VORONOI_VALUE_BORDER = 3,
  VORONOI_VALUE_SECOND = 4,
  VORONOI_VALUE_SPROD = 5,
  VORONOI_VALUE_LAST
} VoronoiValueType;

typedef struct {
  gdouble x, y;
} VoronoiCoords;

typedef struct {
  VoronoiCoords v; /* line equation: v*r == d (XXX: NOT v*r + d == 0) */
  gdouble d;
} VoronoiLine;

typedef struct VoronoiObject {
  VoronoiCoords pos; /* coordinates */
  VoronoiLine rel; /* precomputed coordinates relative to currently processed
                      object and their norm */
  gdouble angle; /* precomputed angle relative to currently processed object
                    (similar as rel) */
  gdouble random;  /* a random number in [0,1], generated to be always the
                      same for the same grid size */
  GSList *ne; /* neighbour list */
} VoronoiObject;

/* convience macro to make the source readable
   use VOBJ(p->next)->angle = M_PI2 */
#define VOBJ(x) ((VoronoiObject*)(x)->data)

#define DOTPROD_SS(a, b) ((a).x*(b).x + (a).y*(b).y)
#define DOTPROD_SP(a, b) ((a).x*(b)->x + (a).y*(b)->y)
#define DOTPROD_PS(a, b) ((a)->x*(b).x + (a)->y*(b).y)
#define DOTPROD_PP(a, b) ((a)->x*(b)->x + (a)->y*(b)->y)

extern gdouble radial_factor;

GSList**        random_grid              (const VoronoiGridType gridtype,
                                          const gint width,
                                          const gint height,
                                          const gdouble side,
                                          const gdouble sigma,
                                          const gdouble lambda,
                                          const gboolean htil,
                                          const gboolean vtil,
                                          const gdouble interstit,
                                          const gdouble vacancies,
                                          gint *wsq,
                                          gint *hsq);
void            destroy_grid_contents    (GSList **squares,
                                          const gint wsq,
                                          const gint hsq);
void            destroy_neighbourhoods   (GSList **squares,
                                          const gint wsq,
                                          const gint hsq);
void            find_voronoi_neighbours  (GSList **squares,
                                          const gint wsq,
                                          const gint hsq);
void            find_voronoi_neighbours_iter  (GSList **squares,
                                               const gint wsq,
                                               const gint hsq,
                                               const gint iter);
/*
void            remove_objects           (GSList **squares,
                                          const gint wsq,
                                          const gint hsq,
                                          gboolean htil,
                                          gboolean vtil);
                                          */
VoronoiObject*  find_owner               (GSList **squares,
                                          const gint wsq,
                                          const gint hsq,
                                          const VoronoiCoords *const point);
VoronoiObject*  move_along_line          (const VoronoiObject *const owner,
                                          const VoronoiCoords *const start,
                                          const VoronoiCoords *const end,
                                          gint *next_safe);
void            neighbourize             (GSList *ne0,
                                          const VoronoiCoords *const center);
void            compute_segment_angles   (GSList *ne0);

typedef gdouble (*VoronoiValueFunc)(const VoronoiCoords* const,
                                    const VoronoiObject* const);

gdouble         value_flat                (const VoronoiCoords *const point,
                                           const VoronoiObject *const owner);
gdouble         value_radial              (const VoronoiCoords *const point,
                                           const VoronoiObject *const owner);
gdouble         value_segmented           (const VoronoiCoords *const point,
                                           const VoronoiObject *const owner);
gdouble         value_border              (const VoronoiCoords *const point,
                                           const VoronoiObject *const owner);
gdouble         value_second              (const VoronoiCoords *const point,
                                           const VoronoiObject *const owner);
gdouble         value_sprod               (const VoronoiCoords *const point,
                                           const VoronoiObject *const owner);

static const VoronoiValueFunc value_function[] = {
  &value_flat,      /* VORONOI_VALUE_FLAT */
  &value_radial,    /* VORONOI_VALUE_RADIAL */
  &value_segmented, /* VORONOI_VALUE_SEGMENTED */
  &value_border,    /* VORONOI_VALUE_BORDER */
  &value_second,    /* VORONOI_VALUE_SECOND */
  &value_sprod,     /* VORONOI_VALUE_SPROD */
};

/* coords_minus(const VoronoiCoords*const, const VoronoiCoords*const) {{{ */
static inline VoronoiCoords
coords_minus(const VoronoiCoords *const a, const VoronoiCoords *const b)
{
  VoronoiCoords z;

  z.x = a->x - b->x;
  z.y = a->y - b->y;

  return z;
} /* }}} */

/* coords_plus(const VoronoiCoords*const, const VoronoiCoords*const) {{{ */
static inline VoronoiCoords
coords_plus(const VoronoiCoords *const a, const VoronoiCoords *const b)
{
  VoronoiCoords z;

  z.x = a->x + b->x;
  z.y = a->y + b->y;

  return z;
} /* }}} */

#endif /* VORONOI_H */
