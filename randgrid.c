/*
 * @(#) $Id: randgrid.c,v 1.28 2004/05/20 17:53:08 yeti Exp $
 * Generate semi-random two-dimensional grids.
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
#include <math.h>
#include <string.h>
#include <glib.h>
#include "voronoi.h"

#define SQRT3 1.73205080756887729352
/* given two integer indices, return centre of grid element coordinates */
typedef VoronoiCoords (*VoronoiGridFunc)(const VoronoiCoords*,
                                         const gint, const gint);

/* areas of objects of side 1 */
static const gdouble shape_area_factor[] = {
  1.0,  /* VORONOI_GRID_SQUARE */
  1.0,  /* VORONOI_GRID_RHOMB */
  SQRT3/4,  /* VORONOI_GRID_TRIANGLE_H */
  SQRT3/4,  /* VORONOI_GRID_TRIANGLE_V */
  3*SQRT3/2,  /* VORONOI_GRID_HEXAGON_H */
  3*SQRT3/2,  /* VORONOI_GRID_HEXAGON_V */
};

static const VoronoiCoords rect_cell_to_axy_ratio[] = {
  { 1.0, 1.0 },  /* VORONOI_GRID_SQUARE */
  { G_SQRT2, G_SQRT2 },  /* VORONOI_GRID_RHOMB */
  { 1.0, SQRT3 },  /* VORONOI_GRID_TRIANGLE_H */
  { SQRT3, 1.0 },  /* VORONOI_GRID_TRIANGLE_V */
  { 3.0, SQRT3 },  /* VORONOI_GRID_HEXAGON_H */
  { SQRT3, 3.0 },  /* VORONOI_GRID_HEXAGON_V */
};

static const gdouble grid_cell_center_distance[] = {
  1.0,  /* VORONOI_GRID_SQUARE */
  1.0,  /* VORONOI_GRID_RHOMB */
  1.0/SQRT3,  /* VORONOI_GRID_TRIANGLE_H */
  1.0/SQRT3,  /* VORONOI_GRID_TRIANGLE_V */
  SQRT3,  /* VORONOI_GRID_HEXAGON_H */
  SQRT3,  /* VORONOI_GRID_HEXAGON_V */
};

static VoronoiCoords grid_point_square(const VoronoiCoords *a,
                                       const gint jx, const gint jy);
static VoronoiCoords grid_point_rhomb(const VoronoiCoords *a,
                                      const gint jx, const gint jy);
static VoronoiCoords grid_point_triangle_h(const VoronoiCoords *a,
                                           const gint jx, const gint jy);
static VoronoiCoords grid_point_triangle_v(const VoronoiCoords *a,
                                           const gint jx, const gint jy);
static VoronoiCoords grid_point_hexagon_h(const VoronoiCoords *a,
                                          const gint jx, const gint jy);
static VoronoiCoords grid_point_hexagon_v(const VoronoiCoords *a,
                                          const gint jx, const gint jy);

static const VoronoiGridFunc shape_grid_function[] = {
  &grid_point_square, /* VORONOI_GRID_SQUARE */
  &grid_point_rhomb, /* VORONOI_GRID_RHOMB */
  &grid_point_triangle_h, /* VORONOI_GRID_TRIANGLE_H */
  &grid_point_triangle_v, /* VORONOI_GRID_TRIANGLE_V */
  &grid_point_hexagon_h, /* VORONOI_GRID_HEXAGON_H */
  &grid_point_hexagon_v, /* VORONOI_GRID_HEXAGON_V */
};

/*
 * +-------+-------+-
 * |       |       |
 * |   o   |   o   |
 * |       |       |
 * +-------+-------+-
 * |       |       |
 * |   o   |   o   |
 *
 * Cell ratio: 1:1
 * Cell area: 1
 * Rectangular cell ratio: 1:1
 */
static VoronoiCoords
grid_point_square(const VoronoiCoords *a,
                  const gint jx, const gint jy)
{
  VoronoiCoords z;

  z.x = (jx + 0.5) * a->x;
  z.y = (jy + 0.5) * a->y;

  return z;
}

static VoronoiCoords
grid_point_rhomb(const VoronoiCoords *a,
                 const gint jx, const gint jy)
{
  VoronoiCoords z;

  z.x = (jx + 0.5*(jy&1) + 0.25)*G_SQRT2 * a->x;
  z.y = (jy + 0.5)/G_SQRT2 * a->y;

  return z;
}

/*
 * +-------+-------+-
 *  \  o  / \  o  / \
 *   \   /   \   /
 * o  \ /  o  \ /  o
 * ----+-------+-----
 * o  / \  o  / \  o
 *   /   \   /   \
 *
 * Cell ratio: 1:sqrt(3)/2
 * Cell area: sqrt(3)/4
 * Rectangular cell ratio: 1:sqrt(3)
 */
static VoronoiCoords
grid_point_triangle_h(const VoronoiCoords *a,
                      const gint jx, const gint jy)
{
  VoronoiCoords z;

  z.x = (jx/2.0 - 0.25) * a->x;
  z.y = ((jy + (jy + (jx&1))/2) + 0.75 - 0.5*(jx&1))/SQRT3 * a->y;

  return z;
}

static VoronoiCoords
grid_point_triangle_v(const VoronoiCoords *a,
                      const gint jx, const gint jy)
{
  VoronoiCoords z;

  z.x = ((jx + (jx + (jy&1))/2) + 0.75 - 0.5*(jy&1))/SQRT3 * a->x;
  z.y = (jy/2.0 - 0.25) * a->y;

  return z;
}

/*
 * +---+    o    +---
 *      \       /
 *       \     /
 *   o    +---+    o
 *       /     \
 *      /       \
 * +---+    o    +---
 *      \       /
 *
 * Cell ratio: 2:sqrt(3)
 * Cell area: 3/2*sqrt(3)
 * Rectangular cell ratio: 3:sqrt(3)
 */
static VoronoiCoords
grid_point_hexagon_h(const VoronoiCoords *a,
                     const gint jx, const gint jy)
{
  VoronoiCoords z;

  z.x = (3*jx - 1.25 + 1.5*(jy&1)) * a->x;
  z.y = SQRT3/2.0*(jy - 0.5) * a->y;

  return z;
}

static VoronoiCoords
grid_point_hexagon_v(const VoronoiCoords *a,
                     const gint jx, const gint jy)
{
  VoronoiCoords z;

  z.x = SQRT3/2.0*(jx - 0.5) * a->x;
  z.y = (3*jy - 1.25 + 1.5*(jx&1)) * a->y;

  return z;
}

static GSList*
objlist_shifted_copy(GSList *ne0, const VoronoiCoords *const shift)
{
  GSList *ne, *ne1;
  VoronoiObject *obj;

  ne = g_slist_copy(ne0);
  for (ne1 = ne; ne1 != NULL; ne1 = ne1->next) {
    obj = g_new(VoronoiObject, 1);
    memcpy(obj, ne0->data, sizeof(VoronoiObject));
    obj->pos = coords_plus(&obj->pos, shift);
    ne1->data = obj;

    ne0 = ne0->next;
  }

  return ne;
}

static gboolean
empty_squares(GSList **squares,
              const gint wsq,
              const gint hsq)
{
  gint jx, jy = 0;
  gint xwsq, xhsq;

  xwsq = wsq + 2*SQBORDER;
  xhsq = hsq + 2*SQBORDER;

  /* check whether all inner squares contain at least one object */
  for (jx = 0; jx < wsq; jx++) {
    for (jy = 0; jy < hsq; jy++) {
      if (squares[xwsq*(jy + SQBORDER) + jx + SQBORDER] == NULL)
        break;
    }
  }

  if (jx == wsq && jy == hsq)
    return FALSE; /* no empty squares (== OK) */

  destroy_grid_contents(squares, wsq, hsq);

  return TRUE; /* some empty squares (== BAD) */
}

static inline void
free_objlist(GSList *l)
{
  GSList *ll;

  for (ll = l; ll; ll = g_slist_next(ll))
    g_free(ll->data);
  g_slist_free(l);
}

static void
random_squarized_points(GSList **squares,
                        const gint wsq,
                        const gint hsq,
                        gint n)
{
  gint still_empty;
  gint xwsq, xhsq;
  gint jx, jy, i, ii;
  gdouble tx, ty;
  VoronoiObject *obj;

  g_assert(n > wsq*hsq);
  xwsq = wsq + 2*SQBORDER;
  xhsq = hsq + 2*SQBORDER;
  tx = (gdouble)xwsq/wsq;
  ty = (gdouble)xhsq/hsq;
  n = tx*ty*n;  /* correct number of points for border */

  /* create points randomly as long as we can,
   * create them carelessly over all extended grid, we will eventually force
   * periodicity later */
  still_empty = xwsq*xhsq;
  for (i = 0; i < n - still_empty; i++) {
    obj = g_new(VoronoiObject, 1);
    obj->pos.x = xwsq*g_rand_double(rng) - SQBORDER;
    obj->pos.y = xhsq*g_rand_double(rng) - SQBORDER;
    obj->random = g_rand_double(rng);
    obj->ne = NULL;

    jx = floor(obj->pos.x) + SQBORDER;
    jy = floor(obj->pos.y) + SQBORDER;

    if (squares[xwsq*jy + jx] == NULL)
      still_empty--;
    squares[xwsq*jy + jx] = g_slist_prepend(squares[xwsq*jy + jx], obj);
  }

  /* the remaining points must be placed into still-empty squares */
  /*fprintf(stderr, "still_empty = %d\n", still_empty);*/
  for (jy = 0; jy < xhsq; jy++) {
    for (jx = 0; jx < xwsq; jx++) {
      ii = jy*xwsq + jx;
      if (squares[ii] == NULL) {
        obj = g_new(VoronoiObject, 1);
        obj->pos.x = g_rand_double(rng) + jx - SQBORDER;
        obj->pos.y = g_rand_double(rng) + jy - SQBORDER;
        obj->random = g_rand_double(rng);
        obj->ne = NULL;
        squares[ii] = g_slist_prepend(squares[ii], obj);
      }
    }
  }

  g_assert(!empty_squares(squares, wsq, hsq));
}

static void
make_grid_tilable(GSList **squares,
                  gint wsq,
                  gint hsq,
                  gboolean htil,
                  gboolean vtil)
{
  VoronoiCoords shift;
  gint jx, jy;
  gint xwsq, xhsq;

  xwsq = wsq + 2*SQBORDER;
  xhsq = hsq + 2*SQBORDER;

  /* finally copy points to the border area, making the grid periodic in
   * requested direction */
  if (vtil) {
    shift.x = 0;
    shift.y = hsq;
    for (jy = 0; jy < SQBORDER; jy++) {
      for (jx = 0; jx < xwsq; jx++) {
        free_objlist(squares[(hsq + jy + SQBORDER)*xwsq + jx]);
        squares[(hsq + jy + SQBORDER)*xwsq + jx]
          = objlist_shifted_copy(squares[(jy + SQBORDER)*xwsq + jx], &shift);
      }
    }
    shift.y = -hsq;
    for (jy = 0; jy < SQBORDER; jy++) {
      for (jx = 0; jx < xwsq; jx++) {
        free_objlist(squares[jy*xwsq + jx]);
        squares[jy*xwsq + jx]
          = objlist_shifted_copy(squares[(jy + hsq)*xwsq + jx], &shift);
      }
    }
  }

  if (htil) {
    shift.x = wsq;
    shift.y = 0;
    for (jy = 0; jy < xhsq; jy++) {
      for (jx = 0; jx < SQBORDER; jx++) {
        free_objlist(squares[jy*xwsq + jx + wsq + SQBORDER]);
        squares[jy*xwsq + jx + wsq + SQBORDER]
          = objlist_shifted_copy(squares[jy*xwsq + jx + SQBORDER], &shift);
      }
    }
    shift.x = -wsq;
    for (jy = 0; jy < xhsq; jy++) {
      for (jx = 0; jx < SQBORDER; jx++) {
        free_objlist(squares[jy*xwsq + jx]);
        squares[jy*xwsq + jx]
          = objlist_shifted_copy(squares[jy*xwsq + jx + wsq], &shift);
      }
    }
  }
}

static void
compute_deformation_map(gdouble *map,
                        gint width,
                        gint height,
                        gdouble lambda)
{
  gint i, j, k, kr;
  gdouble *kernel, *map2;
  gdouble q, ksum;

  /* uncorrelated random numbers */
  for (i = 0; i < width*height; i++)
    map[i] = g_rand_double_range(rng, -1.0, 1.0);

  /* kernel */
  kr = ceil(1.8*lambda);

  /* always make the map cyclic to prevent problems for too large lambda */
  kernel = g_new(gdouble, 2*kr+1);
  for (k = 0; k <= kr; k++)
    kernel[kr + k] = kernel[kr - k] = exp(-2.0*k*k/lambda/lambda);
  kernel += kr;

  ksum = 0.0;
  for (k = -kr; k <= kr; k++) {
    ksum += kernel[k];
    /*fprintf(stderr, "k[%d] = %f\n", k, kernel[k]);*/
  }

  map2 = g_new(gdouble, width*height);

  /* horizontal blur */
  for (i = 0; i < height; i++) {
    gdouble *row = map + i*width;
    gdouble *row2 = map2 + i*width;

    for (j = 0; j < width; j++) {
      gdouble s = 0.0;

      for (k = -kr; k <= kr; k++)
        s += kernel[k]*row[(j + k + 4*width) % width];
      row2[j] = s;
    }
  }

  q = 0.0;
  for (i = 0; i < width*height; i++)
    q += fabs(map[i]);

  /* vertical blur */
  for (j = 0; j < width; j++) {
    gdouble *col = map + j;
    gdouble *col2 = map2 + j;

    for (i = 0; i < height; i++) {
      gdouble s = 0.0;

      for (k = -kr; k <= kr; k++)
        s += kernel[k]*col2[((i + k + 4*height) % height)*width];
      col[i*width] = s;
    }
  }

  q = 0.0;
  for (i = 0; i < width*height; i++)
    q += fabs(map[i]);

  g_free(map2);
  /* normalize */
  q = lambda/(ksum*ksum);
  for (i = 0; i < width*height; i++)
    map[i] *= q;

  q = 0.0;
  for (i = 0; i < width*height; i++)
    q += fabs(map[i]);

  kernel -= kr;
  g_free(kernel);
}

static void
apply_deformation_map(GSList *grid,
                      gint wsq,
                      gint hsq,
                      const gdouble *defmapx,
                      const gdouble *defmapy,
                      VoronoiCoords sigma,
                      gint wmap,
                      gint hmap)
{
  gint jx, jy, jxp, jyp;
  gdouble rx, ry;
  gdouble def;

  while (grid) {
    VoronoiObject *obj = VOBJ(grid);

    jx = floor(obj->pos.x/wsq*wmap);
    jy = floor(obj->pos.y/hsq*hmap);
    rx = obj->pos.x/wsq*wmap - jx;
    ry = obj->pos.y/hsq*hmap - jy;
    jx = (jx + wmap) % wmap;
    jy = (jy + hmap) % hmap;
    jxp = (jx + 1) % wmap;
    jyp = (jy + 1) % hmap;

    def = defmapx[jy*wmap + jx]*(1.0 - rx)*(1.0 - ry)
          + defmapx[jy*wmap + jxp]*rx*(1.0 - ry)
          + defmapx[jyp*wmap + jx]*(1.0 - rx)*ry
          + defmapx[jyp*wmap + jxp]*rx*ry;
    obj->pos.x += sigma.x*def;
    def = defmapy[jy*wmap + jx]*(1.0 - rx)*(1.0 - ry)
          + defmapy[jy*wmap + jxp]*rx*(1.0 - ry)
          + defmapy[jyp*wmap + jx]*(1.0 - rx)*ry
          + defmapy[jyp*wmap + jxp]*rx*ry;
    obj->pos.y += sigma.y*def;

    grid = g_slist_next(grid);
  }
}

/* create almost-prefectly regular grid, deformations are added later */
static GSList*
create_regular_grid(const VoronoiGridType gridtype,
                    const VoronoiCoords *a,
                    const gint wsq,
                    const gint hsq)
{
  GSList *grid;
  gint jx, jy;
  gint xwsq, xhsq;
  VoronoiObject *obj;
  VoronoiCoords z, incr;
  VoronoiGridFunc center_point;

  g_assert(gridtype < G_N_ELEMENTS(shape_grid_function));
  center_point = shape_grid_function[gridtype];

  xwsq = wsq + 2*SQBORDER;
  xhsq = hsq + 2*SQBORDER;
  incr.x = a->x*rect_cell_to_axy_ratio[gridtype].x;
  incr.x = ceil(SQBORDER/incr.x)*incr.x;
  incr.y = a->y*rect_cell_to_axy_ratio[gridtype].y;
  incr.y = ceil(SQBORDER/incr.y)*incr.y;

  grid = NULL;
  for (jy = 0; ; jy++) {
    for (jx = 0; ; jx++) {
      z = center_point(a, jx, jy);
      if (z.x > xwsq)
        break;
      z = coords_minus(&z, &incr);
      /* FIXME add some minimal noise to prevent artifacts and
       * numeric instability */
      z.x += g_rand_double_range(rng, -0.001, 0.001);
      z.y += g_rand_double_range(rng, -0.001, 0.001);

      obj = g_new(VoronoiObject, 1);
      obj->pos = z;
      obj->random = g_rand_double(rng);
      obj->ne = NULL;
      grid = g_slist_prepend(grid, obj);
    }
    if (z.y > xhsq)
      break;
  }

  return grid;
}

static void
add_dislocations(GSList **squares,
                 const gint wsq,
                 const gint hsq,
                 const gint interstit,
                 const gint vacancies)
{
  VoronoiObject *obj;
  GSList *l;
  gint *must_fill;
  gint i, j, jx, jy, k, len;
  gint can_empty, ntofill, misses;
  gint xwsq, xhsq;

  /*fprintf(stderr, "interstit = %d, vacancies = %d\n", interstit, vacancies);*/
  if (!interstit && !vacancies)
    return;

  xwsq = wsq + 2*SQBORDER;
  xhsq = hsq + 2*SQBORDER;

  can_empty = MIN(interstit, vacancies);
  must_fill = g_new(gint, can_empty);
  ntofill = misses = 0;

  /* remove random objects.  we can make some squares empty, however we have
   * to remember to fill them with interstitials later
   * don't try too hard to really create the exact number of vacancies, who
   * knows if it's possible, just try `hard enough'
   * in practice, the !len and len == 1 code paths should be taken *very*
   * rarely */
  i = 0;
  while (i < vacancies && misses <= vacancies) {
    j = g_rand_int_range(rng, 0, xwsq*xhsq);
    len = g_slist_length(squares[j]);

    if (!len) {
      /*fprintf(stderr, "Hit empty (%d,%d)\n", j%xwsq, j/xwsq);*/
      misses++;
      continue;
    }

    if (len == 1) {
      if (ntofill == can_empty) {
        /*fprintf(stderr, "Cannot empty (%d,%d) any more\n", j%xwsq, j/xwsq);*/
        misses++;
        continue;
      }
      /*fprintf(stderr, "Emptying (%d,%d)\n", j%xwsq, j/xwsq);*/
      must_fill[ntofill++] = j;
      g_free(squares[j]->data);
      g_slist_free_1(squares[j]);
      squares[j] = NULL;
      i++;
      continue;
    }

    k = g_rand_int_range(rng, 0, len);
    /*fprintf(stderr, "Removing %d from (%d,%d)\n", k, j%xwsq, j/xwsq);*/
    l = g_slist_nth(squares[j], k);
    g_free(l->data);
    squares[j] = g_slist_delete_link(squares[j], l);
    i++;
  }

  /* add interstitials, filling the empty squares first */
  for (i = 0; i < ntofill; i++) {
    j = must_fill[i];
    jx = j % xwsq;
    jy = j/xwsq;

    obj = g_new(VoronoiObject, 1);
    obj->pos.x = g_rand_double(rng) + jx - SQBORDER;
    obj->pos.y = g_rand_double(rng) + jy - SQBORDER;
    obj->random = g_rand_double(rng);
    obj->ne = NULL;
    squares[j] = g_slist_prepend(squares[j], obj);
  }
  /* the remaining ones are really random */
  for (i = 0; i < interstit - ntofill; i++) {
    obj = g_new(VoronoiObject, 1);
    obj->pos.x = xwsq*g_rand_double(rng) - SQBORDER;
    obj->pos.y = xhsq*g_rand_double(rng) - SQBORDER;
    obj->random = g_rand_double(rng);
    obj->ne = NULL;

    jx = floor(obj->pos.x) + SQBORDER;
    jy = floor(obj->pos.y) + SQBORDER;
    j = xwsq*jy + jx;
    squares[j] = g_slist_prepend(squares[j], obj);
  }
}

static void
squarize_grid(GSList *grid,
              GSList **squares,
              const gint wsq,
              const gint hsq)
{
  GSList *this, *next;
  gint jx, jy;
  gint xwsq, xhsq;
  gdouble t;

  xwsq = wsq + 2*SQBORDER;
  xhsq = hsq + 2*SQBORDER;

  for (this = grid; this != NULL; this = next) {
    VoronoiObject *obj = VOBJ(this);

    jx = t = floor(obj->pos.x) + SQBORDER;
    if (t == obj->pos.x + SQBORDER)
      obj->pos.x += EPS; /* degenerate case */

    jy = t = floor(obj->pos.y) + SQBORDER;
    if (t == obj->pos.y + SQBORDER)
      obj->pos.y += EPS; /* degenerate case */

    next = this->next;
    if (jx >= 0 && jx < xwsq && jy >= 0 && jy < xhsq) {
      /* relink it from point grid to the appropriate square */
      this->next = squares[xwsq*jy + jx];
      squares[xwsq*jy + jx] = this;
    }
    else {
      g_free(this->data);
      g_slist_free_1(this); /* discard it */
    }
  }
}

/* does NOT destroy neighbourhoods */
void
destroy_grid_contents(GSList **squares,
                      const gint wsq,
                      const gint hsq)
{
  gint xwsq, xhsq, i;

  xwsq = wsq + 2*SQBORDER;
  xhsq = hsq + 2*SQBORDER;
  for (i = 0; i < xwsq*xhsq; i++)
    free_objlist(squares[i]);
}

/**
 * random_grid:
 * @gridtype: Grid type to construct.
 * @width: Grid width in some units (e.g., pixels).
 * @height: Grid height in some units (e.g., pixels).
 * @side: Grid element side in the same units as @width, @height.
 * @sigma: Random deformation sigma.  Dimensionless, relative to @side.
 * @lambda: Autocorrelation length of random deformation.  Dimensionless,
 *          relative to @side.
 * @htil: Whether it should be horizontally tileable.
 * @vtil: Whether it should be vertically tileable.
 * @wsq: Grid width (in squares, without border) will be stored here.
 * @hsq: Grid height (in squares, without border) will be stored here.
 *
 * Create a random grid if type @gridtype.
 *
 * The gird is already squarized, i.e., it is transformed into a list of
 * squares containing corresponding grid points.  The size of each square
 * is 1x1, @wsq/@hsq ratio is approximately equal to @width/@height.
 *
 * Returns: The grid as an array of @wsq+2*%SQBORDER times @hsq+2*%SQBORDER
 *          #VoronoiObject lists.
 **/
GSList**
random_grid(const VoronoiGridType gridtype,
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
            gint *hsq)
{
  GSList *grid;
  GSList **squares;
  VoronoiCoords ncell, a, defsigma;
  gdouble slambda;
  gdouble *defmapx, *defmapy;
  gint xwsq, xhsq;
  gint iter, n;
  gint wmap, hmap;

  g_return_val_if_fail(gridtype == VORONOI_GRID_RANDOM
                       || gridtype <= G_N_ELEMENTS(shape_grid_function),
                       NULL);
  /* compute square size trying to keep the aspect ratio */
  if (width < height) {
    *wsq = width/(sqrt(7.0)*side);
    a.x = width/(gdouble)*wsq;
    *hsq = floor(height/a.x + 0.5);
  }
  else {
    *hsq = height/(sqrt(7.0)*side);
    a.y = height/(gdouble)*hsq;
    *wsq = floor(width/a.y + 0.5);
  }
  n = ceil((width*height)/(side*side));
  /*fprintf(stderr, "n = %d\n", n);*/
  xwsq = *wsq + 2*SQBORDER;
  xhsq = *hsq + 2*SQBORDER;

  squares = g_new0(GSList*, xwsq*xhsq);

  if (gridtype == VORONOI_GRID_RANDOM) {
    /* assure different cell sizes result in different random grids, even if
     * hsq and wsq are the same, it only fixes preview user experience,
     * nothing more */
    for (iter = 0; iter < (gint)(83*side) % 37; iter++)
      g_rand_int(rng);

    random_squarized_points(squares, *wsq, *hsq, n);
    make_grid_tilable(squares, *wsq, *hsq, htil, vtil);
    radial_factor = sqrt((width*height)/(gdouble)((*wsq)*(*hsq)))/side;
    return squares;
  }

  /* compute cell sizes, must take tileability into account
   * ncell is the number of *rectangular* cells */
  ncell.x = width/side*sqrt(shape_area_factor[gridtype])
            /rect_cell_to_axy_ratio[gridtype].x;
  ncell.y = height/side*sqrt(shape_area_factor[gridtype])
            /rect_cell_to_axy_ratio[gridtype].y;
  if (htil && vtil) {
    /* tilable in both directions => deformation is inevitable */
    if (ncell.x < ncell.y) {
      /* use x to compute cell width, deform cell height */
      a.y = a.x = width/ceil(ncell.x)/rect_cell_to_axy_ratio[gridtype].x;
      ncell.y = height/a.y/rect_cell_to_axy_ratio[gridtype].y;
      a.y = height/floor(ncell.y + 0.5)/rect_cell_to_axy_ratio[gridtype].y;
    }
    else {
      /* use y to compute cell height, deform cell width */
      a.x = a.y = height/ceil(ncell.y)/rect_cell_to_axy_ratio[gridtype].y;
      ncell.x = width/a.x/rect_cell_to_axy_ratio[gridtype].x;
      a.x = width/floor(ncell.x + 0.5)/rect_cell_to_axy_ratio[gridtype].x;
    }
  }
  else if (htil) {
    /* just htil => compute vertical cell size to match horizontal */
      a.y = a.x = width/ceil(ncell.x)/rect_cell_to_axy_ratio[gridtype].x;
  }
  else if (vtil) {
    /* just vtil => compute horizontal cell size to match vertical */
      a.x = a.y = height/ceil(ncell.y)/rect_cell_to_axy_ratio[gridtype].y;
  }
  else {
    /* no tilability */
      a.x = width/ncell.x/rect_cell_to_axy_ratio[gridtype].x;
      a.y = height/ncell.y/rect_cell_to_axy_ratio[gridtype].y;
  }
  /*fprintf(stderr, "a = (%f, %f), %f\n", a.x, a.y, a.x/a.y);*/

  /* compensate for square deformation */
  a.x *= *wsq/(gdouble)width;
  a.y *= *hsq/(gdouble)height;
  radial_factor = 1.0/sqrt(a.x*a.y*grid_cell_center_distance[gridtype]);

  /*fprintf(stderr, "lambda = %f\n", lambda);*/
  wmap = MAX(1.2*rect_cell_to_axy_ratio[gridtype].x*ncell.x, 6);
  hmap = MAX(1.2*rect_cell_to_axy_ratio[gridtype].y*ncell.y, 6);
  defmapx = g_new(gdouble, wmap*hmap);
  defmapy = g_new(gdouble, wmap*hmap);
  slambda = lambda*sqrt(wmap*hmap/(gdouble)n);
  /*
  fprintf(stderr, "ncell = (%f, %f)\n", ncell.x, ncell.y);
  fprintf(stderr, "a = (%f, %f)\n", a.x, a.y);
  fprintf(stderr, "wmap = %d, hmap = %d, slambda = %f\n", wmap, hmap, slambda);
  */

  defsigma.x = a.x*sigma*grid_cell_center_distance[gridtype];
  defsigma.y = a.y*sigma*grid_cell_center_distance[gridtype];
  /*fprintf(stderr, "sigma = (%f, %f)\n", defsigma.x, defsigma.y);*/
  n = (xwsq*xhsq)/(a.x*a.y*shape_area_factor[gridtype]);
  iter = 0;
  do {
    grid = create_regular_grid(gridtype, &a, *wsq, *hsq);
    compute_deformation_map(defmapx, wmap, hmap, slambda);
    compute_deformation_map(defmapy, wmap, hmap, slambda);
    apply_deformation_map(grid, *wsq, *hsq,
                          defmapx, defmapy, defsigma, wmap, hmap);
    squarize_grid(grid, squares, *wsq, *hsq);
    add_dislocations(squares, *wsq, *hsq, n*interstit, n*vacancies);
    make_grid_tilable(squares, *wsq, *hsq, htil, vtil);
    if (iter++ == 20) {
      g_error("Cannot create grid. If you can reproduce this, report bug "
              "to <yeti@physics.muni.cz>, please include grid generator "
              "settings.");
    }
  } while (empty_squares(squares, *wsq, *hsq));

  g_free(defmapx);
  g_free(defmapy);

  return squares;
}
