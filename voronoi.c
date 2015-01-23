/*
 * @(#) $Id: voronoi.c,v 1.39 2004/05/20 18:03:17 yeti Exp $
 * This is a plugin for the GIMP.
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

/*
 * This plug-in generates 2D Voronoi diagrams of various semi-random sets of
 * points and knows some nice ways how to color them.
 *
 * I cheat, because it's not possible to generate a Voronoi diagram in O(n)
 * time fairly.  However your chance to catch me is less than 1:1000
 * in the worst case, and negligible normally.
 *
 * Please send patches, suggestions, and other invectives to me:
 * <yeti@physics.muni.cz>.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#else
   /* To include libgimp/stdplugins-intl.h */
#  define GETTEXT_PACKAGE "gimp20"
#endif /* HAVE_CONFIG_H */

#include <stdlib.h>
#include <string.h>

#include <gtk/gtk.h>

#include <libgimp/gimp.h>
#include <libgimp/gimpui.h>
#include "libgimp/stdplugins-intl.h"

#include "voronoi.h"

#define PROCEDURE_NAME   "voronoi"
#define DATA_KEY_VALS    "plug_in_voronoi"

/* we are not 100% 16-bit-ready, but partially yes (almost everything is done
   on floats anyway) */
#define CHANNEL_MAX_VALUE     255
#define BITS_PER_SAMPLE         8

#define PREVIEW_SIZE          128
#define CELL_SIZE_WIDTH        84

/* Don't change, though it's ugly */
#define SCALE_WIDTH 201

/* XXX: we cannot store true strings as parameters, hope this is enough */
#define GRADIENT_NAME_SIZE    256

/* generally, we need much more than 256
   note we do NOT do our own subsampling, just rounding */
#define GRADIENT_SAMPLE_SIZE  8192

#define PREVIEW_CELL_SIZE_MAX   48

enum {
  VORONOI_RESPONSE_REVERT = 1,
  VORONOI_RESPONSE_RESET,
  VORONOI_RESPONSE_ABOUT
};

enum {
  VORONOI_RECOMPUTE_GRID =  1 << 0,
  VORONOI_RECOMPUTE_COLOR = 1 << 1
};

typedef enum {
  VORONOI_CM_GRAY = 0,
  VORONOI_CM_GRADIENT = 1,
} VoronoiColoringMethod;

typedef struct {
  gint                   seed;
  gint                   celltype;
  gdouble                cellsize;
  gdouble                sigma;
  gdouble                lambda;
  gboolean               htil;
  gboolean               vtil;
  gdouble                valscale;
  gdouble                valshift;
  gdouble                interstit;
  gdouble                vacancies;
  gdouble                weights[VORONOI_VALUE_LAST];
  VoronoiColoringMethod  cm;
  gchar                  grad[GRADIENT_NAME_SIZE];
  gboolean               grad_invert;
  gboolean               random_seed;
} VoronoiValues;

/* local prototypes */
static void query(void);
static void run                                (const gchar *name,
                                                gint nparams,
                                                const GimpParam *param,
                                                gint *nreturn_vals,
                                                GimpParam  **return_vals);
static gint voronoi_dialog                     (void);
static void voronoi_refresh_controls           (const VoronoiValues *new_vvals);
static void voronoi_celltype_changed           (GtkWidget *widget,
                                                gpointer data);
static void voronoi_celltype_update_dependants (void);
static void voronoi_cm_changed                 (GtkWidget *widget,
                                                gpointer data);
static void voronoi_gradient_changed           (const gchar *name,
                                                gint width,
                                                const gdouble *grad_data,
                                                gboolean dialog_closing,
                                                gpointer user_data);
static void voronoi_cellsize_changed           (GtkAdjustment *adj);
static void voronoi                            (GimpDrawable *drawable,
                                                guint recompute);
static void init_voronoi                       (GimpDrawable *drawable);
static void commit                             (GimpDrawable *drawable,
                                                gint wsq,
                                                gint hsq,
                                                const gint chmax);
static void commit_one_tile                    (GimpPixelRgn *rgn,
                                                gint wsq,
                                                gint hsq,
                                                const gint chmax,
                                                const gint gradwidth);
static void end_voronoi                        (GimpDrawable *drawable);
static void about_dialog                       (GtkWidget *parent);

GRand *rng = NULL;  /* global random generator, used in randgrid */

/***** Local vars *****/

typedef struct {
  GtkWidget *seed;
  GtkWidget *celltype;
  GtkObject *cellsize;
  GtkObject *sigma;
  GtkObject *lambda;
  GtkWidget *htil;
  GtkWidget *vtil;
  GtkObject *valscale;
  GtkObject *valshift;
  GtkObject *interstit;
  GtkObject *vacancies;
  GtkObject *weights[VORONOI_VALUE_LAST];
  GtkWidget *cm;
  GtkWidget *grad;
  GtkWidget *grad_invert;
  GtkWidget *grad_label;
  GtkWidget *preview_frame;
  GtkWidget *deform_frame;
} VoronoiControls;

static VoronoiControls vctrl;

GimpPlugInInfo PLUG_IN_INFO = {
  NULL,  /* init_proc  */
  NULL,  /* quit_proc  */
  query, /* query_proc */
  run,   /* run_proc   */
};

static const VoronoiValues vvals_defaults = {
  11,    /* seed */
  VORONOI_GRID_SQUARE,     /* celltype */
  20,     /* cellsize */
  0.15,  /* sigma */
  0.7,  /* lambda */
  FALSE,  /* htil */
  FALSE,  /* vtil */
  0.0,  /* valscale */
  0.0,  /* valshift */
  0.0,  /* interstit */
  0.0,  /* vacancies */
  { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, },
  VORONOI_CM_GRAY,   /* cm */
  "",   /* grad */
  FALSE,  /* grad_invert */
  FALSE,  /* random_seed */
};

static const gchar *value_description[VORONOI_VALUE_LAST] = {
  "R_andom solid color:",
  "_Radial distance:",
  "Se_gment distance:",
  "_Border distance:",
  "_Second nearest distance:",
  "Scalar _product:",
};

static VoronoiValues vvals; /* = vvals_defaults; */
static VoronoiValues vvals_old; /* for revert */

static GtkWidget *preview_image = NULL;
static GdkPixbuf *preview = NULL;
static gboolean preview_mode;
static gboolean is_rgb;
static GimpDrawable *drawable;
static gdouble *gradient = NULL;
static gint gradient_p_width = -1;
static gboolean in_param_reset = TRUE;

/*
 * Some globals to save passing too many paramaters that don't change.
 */

static gint     ix1, iy1, ix2, iy2;     /* Selected image size. */
static gint     bpp, has_alpha, alpha;
static gint     voronoi_width, voronoi_height, wsq, hsq;
static glong    max_progress, progress;
static gdouble  value_shift;
static gdouble  value_scale;
static GSList **squares = NULL;

/***** Functions *****/

MAIN()

static void
query(void)
{
  static GimpParamDef args[] = {
    { GIMP_PDB_INT32, "run_mode", "Interactive, non-interactive" },
    { GIMP_PDB_IMAGE, "image", "Input image (unused)" },
    { GIMP_PDB_DRAWABLE, "drawable", "Input drawable" },
    { GIMP_PDB_INT32, "seed", "Random seed" },
    { GIMP_PDB_INT32, "celltype", "Cell type "
                                  "{ RANDOM (-1), SQUARE (0), "
                                  "RHOMB (1), TRIANGLE_H (2), TRIANGLE_V (3) "
                                  "HEXAGON_H (4), HEXAGON_V (5) }" },
    { GIMP_PDB_FLOAT, "cellsize", "Average cell size" },
    { GIMP_PDB_FLOAT, "sigma", "Grid randomness "
                               "(unused for RANDOM)" },
    { GIMP_PDB_FLOAT, "lambda", "Autocorrelation length "
                                "(unused for RANDOM)" },
    { GIMP_PDB_INT32, "htil", "Horizontal tileability" },
    { GIMP_PDB_INT32, "vtil", "Vertical tileability" },
    { GIMP_PDB_FLOAT, "valscale", "Value rescale (log)" },
    { GIMP_PDB_FLOAT, "valshift", "Value shift" },
    { GIMP_PDB_INT32, "cm", "Coloring method: "
                            "{ INDEPENDENT_RGB (0), GRADIENT (1) }" },
    { GIMP_PDB_FLOAT, "interstit", "Interstitials percentage "
                                   "(unused for RANDOM)" },
    { GIMP_PDB_FLOAT, "vacancies", "Vacancies percentage "
                                   "(unused for RANDOM)" },
    { GIMP_PDB_FLOAT, "w_flat", "Random solid color weight" },
    { GIMP_PDB_FLOAT, "w_radial", "Radial distace weight" },
    { GIMP_PDB_FLOAT, "w_segment", "Segment distace weight" },
    { GIMP_PDB_FLOAT, "w_border", "Border distace weight" },
    { GIMP_PDB_FLOAT, "w_second", "Second nearest distace weight" },
    { GIMP_PDB_FLOAT, "w_sprod", "Scalar product weight" },
    { GIMP_PDB_INT32, "grad_invert", "Invert gradient (or gray scale)? " },
    { GIMP_PDB_STRING, "grad", "Gradient used for coloring "
                               "(when cm is GRADIENT)" },
  };

  gimp_install_procedure
    (DATA_KEY_VALS,
     "Render a Voronoi diagram to specified drawable",
     "This plugin creates Voronoi diagram of semi-random points, "
     "and colors it somewhow.",
     "Yeti",
     "Yeti",
     "Oct 2002",
     N_("<Image>/Filters/Render/Pattern/Voronoi..."),
     "RGB*, GRAY*",
     GIMP_PLUGIN,
     G_N_ELEMENTS(args), 0, args, NULL);
}

static void
run(const gchar *name,
    gint nparams,
    const GimpParam *param,
    gint *nreturn_vals,
    GimpParam **return_vals)
{
  static GimpParam values[1];
  GimpRunMode run_mode;
  GimpPDBStatusType status = GIMP_PDB_SUCCESS;
  guint i;

  run_mode = param[0].data.d_int32;

  /*  Initialize i18n support  */
  bindtextdomain(GETTEXT_PACKAGE, LOCALEDIR);
#ifdef HAVE_BIND_TEXTDOMAIN_CODESET
  bind_textdomain_codeset(GETTEXT_PACKAGE, "UTF-8");
#endif
  textdomain(GETTEXT_PACKAGE);

  *nreturn_vals = 1;
  *return_vals = values;

  values[0].type = GIMP_PDB_STATUS;
  values[0].data.d_status = status;

  /*  Get the specified drawable  */
  drawable = gimp_drawable_get(param[2].data.d_drawable);
  is_rgb = gimp_drawable_is_rgb(drawable->drawable_id);

  vvals = vvals_defaults;
  rng = g_rand_new();
  switch (run_mode) {
    case GIMP_RUN_INTERACTIVE:
    /*  Possibly retrieve data  */
    gimp_get_data(DATA_KEY_VALS, &vvals);

    /*  First acquire information with a dialog  */
    vvals_old = vvals;

    if (!voronoi_dialog()) {
      gimp_drawable_detach(drawable);
      return;
    }
    break;

    case GIMP_RUN_NONINTERACTIVE:
    /*  Make sure all the arguments are there!  */
    if (nparams != 15 + VORONOI_VALUE_LAST + 2)
      status = GIMP_PDB_CALLING_ERROR;
    else {
      vvals.seed = (gint)param[3].data.d_int32;
      vvals.celltype = (VoronoiGridType)param[4].data.d_int32;
      vvals.cellsize = (gdouble)param[5].data.d_float;
      vvals.sigma = (gdouble)param[6].data.d_float;
      vvals.lambda = (gdouble)param[7].data.d_float;
      vvals.htil = (gboolean)param[8].data.d_int32;
      vvals.vtil = (gboolean)param[9].data.d_int32;
      vvals.valscale = (gdouble)param[10].data.d_float;
      vvals.valshift = (gdouble)param[11].data.d_float;
      vvals.cm = (VoronoiColoringMethod)param[12].data.d_int32;
      vvals.interstit = (gdouble)param[13].data.d_float;
      vvals.vacancies = (gdouble)param[14].data.d_float;
      for (i = 0; i < VORONOI_VALUE_LAST; i++)
        vvals.weights[i] = (gdouble)param[15 + i].data.d_float;
      i = 15 + VORONOI_VALUE_LAST;
      vvals.grad_invert = (gdouble)param[i].data.d_int32;
      strncpy(vvals.grad, (gchar*)param[i+1].data.d_string,
              GRADIENT_NAME_SIZE-1);
      vvals.grad[GRADIENT_NAME_SIZE-1] = '\0';
    }
    break;

    case GIMP_RUN_WITH_LAST_VALS:
    /*  Possibly retrieve data  */
    gimp_get_data(DATA_KEY_VALS, &vvals);
    break;

    default:
    break;
  }

  if (status == GIMP_PDB_SUCCESS) {
    /*  Make sure that the drawable is gray or RGB color  */
    if (is_rgb || gimp_drawable_is_gray(drawable->drawable_id)) {
      gimp_progress_init(_("Voronoi..."));

      preview_mode = FALSE;
      voronoi(drawable, VORONOI_RECOMPUTE_GRID);

      if (run_mode != GIMP_RUN_NONINTERACTIVE)
        gimp_displays_flush();

      /*  Store data  */
      if (run_mode == GIMP_RUN_INTERACTIVE
          || run_mode == GIMP_RUN_WITH_LAST_VALS)
        gimp_set_data(DATA_KEY_VALS, &vvals, sizeof(VoronoiValues));
    }
    else
      status = GIMP_PDB_EXECUTION_ERROR;
  }

  values[0].data.d_status = status;
  gimp_drawable_detach(drawable);
}

/* remove some source code redundancies {{{ */
static inline void
yeti_table_attach_to_row(GtkWidget *table,
                         GtkWidget *child,
                         const gint row)
{
  gtk_table_attach(GTK_TABLE(table), child, 0, GTK_TABLE(table)->ncols,
                   row, row+1, GTK_EXPAND|GTK_FILL, 0, 0, 0);
}

static inline void
yeti_table_attach_single_cell(GtkWidget *table,
                              GtkWidget *child,
                              const gint column,
                              const gint row)
{
  gtk_table_attach(GTK_TABLE(table), child, column, column+1, row, row+1,
                   GTK_EXPAND|GTK_FILL, 0, 0, 0);
}

static inline GtkWidget*
yeti_label_new_in_table(const gchar *name,
                        GtkWidget *table,
                        const gint row,
                        const gboolean sensitive)
{
  GtkWidget *label;

  label = gtk_label_new_with_mnemonic(name);
  gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
  yeti_table_attach_single_cell(table, label, 0, row);
  gtk_widget_set_sensitive(label, sensitive);

  return label;
}

static inline GtkWidget*
yeti_frame_new_in_box(const gchar *name,
                      GtkWidget *box,
                      const gboolean expand,
                      const gboolean fill)
{
  GtkWidget *frame;

  frame = gtk_frame_new(name);
  gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_ETCHED_IN);
  gtk_box_pack_start(GTK_BOX(box), frame, expand, fill, 0);

  return frame;
}

static inline GtkWidget*
yeti_table_new_in_frame(const gint cols,
                        const gint rows,
                        GtkWidget *frame)
{
  GtkWidget *table;
  GtkWidget *align;

  align = gtk_alignment_new(0.5, 0.0, 1.0, 0.0);
  gtk_container_set_border_width(GTK_CONTAINER(align), 0);
  gtk_container_add(GTK_CONTAINER(frame), align);

  table = gtk_table_new(cols, rows, FALSE);
  gtk_table_set_col_spacings(GTK_TABLE(table), 4);
  gtk_table_set_row_spacings(GTK_TABLE(table), 2);
  gtk_container_set_border_width(GTK_CONTAINER(table), 4);
  gtk_container_add(GTK_CONTAINER(align), table);

  return table;
}

static inline GtkWidget*
yeti_check_button_new_with_label(const gchar *name,
                                 gint *value,
                                 GCallback cb,
                                 gpointer data)
{
  GtkWidget *check;

  check = gtk_check_button_new_with_mnemonic(name);
  g_signal_connect(check, "toggled",
                   G_CALLBACK(gimp_toggle_button_update), value);
  g_signal_connect_swapped(check, "toggled", cb, data);
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(check), *value);

  return check;
}

static inline GtkWidget*
yeti_preview_frame_new_in_box(GtkWidget *hbox,
                              GtkWidget **preview_frame)
{
  GtkWidget *frame;
  GtkWidget *align;
  GtkWidget *image;
  GdkPixbuf *pixbuf;

  frame = yeti_frame_new_in_box(_("Preview"), hbox, FALSE, FALSE);
  *preview_frame = frame;

  align = gtk_alignment_new(0.5, 0.5, 0.0, 0.0);
  gtk_container_set_border_width(GTK_CONTAINER(align), 4);
  gtk_container_add(GTK_CONTAINER(frame), align);

  frame = gtk_frame_new(NULL);
  gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_IN);
  gtk_container_add(GTK_CONTAINER(align), frame);

  pixbuf = gdk_pixbuf_new(GDK_COLORSPACE_RGB,
                          has_alpha,
                          BITS_PER_SAMPLE,
                          PREVIEW_SIZE,
                          PREVIEW_SIZE);
  gdk_pixbuf_fill(pixbuf, 0);

  image = gtk_image_new_from_pixbuf(pixbuf);
  g_object_unref(pixbuf);
  gtk_container_add(GTK_CONTAINER(frame), image);

  return frame;
}

static inline GtkObject*
yeti_scale_entry_new_double(const gchar *name,
                            GtkWidget *table,
                            const gint row,
                            double *value,
                            const double min,
                            const double max,
                            const double step,
                            const double page,
                            const gint digits,
                            GCallback cb,
                            gpointer data)
{
  GtkObject *adj;

  *value = CLAMP(*value, min, max);
  adj = gimp_scale_entry_new(GTK_TABLE(table), 0, row, name, SCALE_WIDTH, 0,
                             *value, min, max, step, page, digits,
                             TRUE, 0, 0, NULL, NULL);
  g_signal_connect(adj, "value_changed",
                   G_CALLBACK(gimp_double_adjustment_update),
                   value);
  g_signal_connect_swapped(adj, "value_changed", cb, data);

  return adj;
}

static inline void
yeti_progress_update(gint p,
                     gint max)
{
  double r;

  if (preview_mode)
    return;

  r = (double)p/max;
  gimp_progress_update(r);
}
/* }}} */

static void
voronoi_g(GimpDrawable *drawable)
{
  voronoi(drawable, VORONOI_RECOMPUTE_GRID);
}

static void
voronoi_c(GimpDrawable *drawable)
{
  voronoi(drawable, VORONOI_RECOMPUTE_COLOR);
}

static gint
voronoi_dialog(void)
{
  GtkWidget *dlg;
  GtkWidget *main_vbox;
  GtkWidget *hbox, *vbox;
  GtkWidget *frame;
  GtkWidget *table;
  GtkWidget *notebook;
  guint i;
  gint run;
  gboolean ok;

  preview_mode = TRUE;
  in_param_reset = TRUE;
  gimp_ui_init(PROCEDURE_NAME, TRUE);

  dlg = gimp_dialog_new(_("Voronoi"), PROCEDURE_NAME,
                        NULL, 0,
                        gimp_standard_help_func, "filters/voronoi.html",

                        _("About"), VORONOI_RESPONSE_ABOUT,
                        _("Revert"), VORONOI_RESPONSE_REVERT,
                        GIMP_STOCK_RESET, VORONOI_RESPONSE_RESET,
                        GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                        GTK_STOCK_OK, GTK_RESPONSE_OK,
                        NULL);

  main_vbox = gtk_vbox_new(FALSE, 4);
  gtk_container_set_border_width(GTK_CONTAINER(main_vbox), 6);
  gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg)->vbox), main_vbox, TRUE, TRUE, 0);

  /* top part with preview and common options */
  hbox = gtk_hbox_new(FALSE, 8);
  gtk_box_pack_start(GTK_BOX(main_vbox), hbox, TRUE, TRUE, 0);

  /* make a nice preview frame */
  frame = yeti_preview_frame_new_in_box(hbox, &vctrl.preview_frame);
  preview_image = GTK_BIN(frame)->child;
  preview = gtk_image_get_pixbuf(GTK_IMAGE(preview_image));

  /* graphic options */
  frame = yeti_frame_new_in_box(_("Graphic Options"), hbox, TRUE, TRUE);
  table = yeti_table_new_in_frame(3, 5, frame);

  /* horizontal tileability */
  vctrl.htil = yeti_check_button_new_with_label(_("_Horizontally tileable"),
                                                &vvals.htil,
                                                G_CALLBACK(voronoi_g),
                                                drawable);
  yeti_table_attach_to_row(table, vctrl.htil, 0);

  /* vertical tileability */
  vctrl.vtil = yeti_check_button_new_with_label(_("_Vertically tileable"),
                                                &vvals.vtil,
                                                G_CALLBACK(voronoi_g),
                                                drawable);
  yeti_table_attach_to_row(table, vctrl.vtil, 1);

  /* gradient, only if rgb */
  if (!is_rgb)
    vvals.cm = VORONOI_CM_GRAY;
  /* synchronize vvals.grad with current gradient */
  if (!strlen(vvals.grad) || !gimp_gradients_set_gradient(vvals.grad)) {
    strncpy(vvals.grad, gimp_gradients_get_gradient(), GRADIENT_NAME_SIZE-1);
    vvals.grad[GRADIENT_NAME_SIZE-1] = '\0';
  }

  vctrl.grad = gimp_gradient_select_widget_new("Voronoi Gradient", vvals.grad,
                                               voronoi_gradient_changed, NULL);
  yeti_table_attach_single_cell(table, vctrl.grad, 1, 3);
  gtk_widget_set_sensitive(vctrl.grad, vvals.cm == VORONOI_CM_GRADIENT);

  vctrl.grad_label = yeti_label_new_in_table(_("_Gradient:"), table, 3,
                                             vvals.cm == VORONOI_CM_GRADIENT);
  gtk_label_set_mnemonic_widget(GTK_LABEL(vctrl.grad_label), vctrl.grad);

  /* coloring function */
  vctrl.cm
    = gimp_option_menu_new2(FALSE, G_CALLBACK(voronoi_cm_changed),
                            vctrl.grad,
                            GUINT_TO_POINTER(vvals.cm),
                            _("Grayscale"),
                            GUINT_TO_POINTER(VORONOI_CM_GRAY), NULL,
                            _("Gradient"),
                            GUINT_TO_POINTER(VORONOI_CM_GRADIENT), NULL,
                            NULL);
  gimp_table_attach_aligned(GTK_TABLE(table), 0, 2,
                            _("Coloring _method:"), 1.0, 0.5,
                            vctrl.cm, 4, TRUE);
  if (!is_rgb)
    gtk_widget_set_sensitive(vctrl.cm, FALSE);

  vctrl.grad_invert = yeti_check_button_new_with_label(_("Invert"),
                                                       &vvals.grad_invert,
                                                       G_CALLBACK(voronoi_c),
                                                       drawable);
  yeti_table_attach_single_cell(table, vctrl.grad_invert, 1, 4);

  notebook = gtk_notebook_new();
  gtk_box_pack_start(GTK_BOX(main_vbox), notebook, TRUE, TRUE, 2);

  /* grid generator */
  vbox = gtk_vbox_new(FALSE, 0);
  gtk_container_set_border_width (GTK_CONTAINER (vbox), 4);
  gtk_notebook_append_page(GTK_NOTEBOOK(notebook), vbox,
                           gtk_label_new(_("Grid Generator")));

  frame = yeti_frame_new_in_box(_("Cells"), vbox, TRUE, TRUE);
  table = yeti_table_new_in_frame(3, 3, frame);

  /* seed/time */
  vctrl.seed = gimp_random_seed_new(&vvals.seed, &vvals.random_seed);
  gimp_table_attach_aligned(GTK_TABLE(table), 0, 0,
                            _("Random Seed:"), 1.0, 0.5,
                            vctrl.seed, 1, TRUE);
  g_signal_connect_swapped(GIMP_RANDOM_SEED_SPINBUTTON_ADJ(vctrl.seed),
                           "value_changed",
                           G_CALLBACK(voronoi_g),
                           drawable);

  /* cell size */
  {
    gint x1, x2, y1, y2;
    gdouble max;

    gimp_drawable_mask_bounds(drawable->drawable_id, &x1, &y1, &x2, &y2);
    max = MIN(x2 - x1, y2 - y1)/sqrt(7.0)/2.01;
    vctrl.cellsize = yeti_scale_entry_new_double(_("Average cell _size"),
                                                 table, 2, &vvals.cellsize,
                                                 3.0, max, 1, 10, 1,
                                                 G_CALLBACK(voronoi_g),
                                                 drawable);
    g_signal_connect(vctrl.cellsize, "value_changed",
                     G_CALLBACK(voronoi_cellsize_changed), &vctrl);
  }

  /* grid (in fact neighbourhood) type */
  vctrl.celltype
    = gimp_option_menu_new2(FALSE, G_CALLBACK(voronoi_celltype_changed), NULL,
                            GUINT_TO_POINTER(vvals.celltype),
                            _("Random"),
                            GUINT_TO_POINTER(VORONOI_GRID_RANDOM), NULL,
                            _("Square"),
                            GUINT_TO_POINTER(VORONOI_GRID_SQUARE), NULL,
                            _("Rhomb"),
                            GUINT_TO_POINTER(VORONOI_GRID_RHOMB), NULL,
                            _("Triangle (horizontal)"),
                            GUINT_TO_POINTER(VORONOI_GRID_TRIANGLE_H), NULL,
                            _("Triangle (vertical)"),
                            GUINT_TO_POINTER(VORONOI_GRID_TRIANGLE_V), NULL,
                            _("Hexagon (horizontal)"),
                            GUINT_TO_POINTER(VORONOI_GRID_HEXAGON_H), NULL,
                            _("Hexagon (vertical)"),
                            GUINT_TO_POINTER(VORONOI_GRID_HEXAGON_V), NULL,
                            NULL);
  gimp_table_attach_aligned(GTK_TABLE(table), 0, 3,
                            _("Cell t_ype:"), 1.0, 0.5,
                            vctrl.celltype, 2, TRUE);

  /* noise */
  frame = yeti_frame_new_in_box(_("Deformation"), vbox, TRUE, TRUE);
  vctrl.deform_frame = frame;
  table = yeti_table_new_in_frame(3, 4, frame);

  /* sigma, lambda */
  vctrl.sigma = yeti_scale_entry_new_double(_("R_andomness:"),
                                            table, 0, &vvals.sigma,
                                            0.0, 2.5, 0.01, 0.1, 3,
                                            G_CALLBACK(voronoi_g),
                                            drawable);
  vctrl.lambda = yeti_scale_entry_new_double(_("Corre_lation:"),
                                               table, 1, &vvals.lambda,
                                               0.5, 8.0, 0.1, 1.0, 3,
                                               G_CALLBACK(voronoi_g),
                                               drawable);
  /* dislocations */
  vctrl.interstit = yeti_scale_entry_new_double(_("Interstitials:"),
                                                table, 2, &vvals.interstit,
                                                0.0, 0.3, 0.001, 0.01, 3,
                                                G_CALLBACK(voronoi_g),
                                                drawable);
  vctrl.vacancies = yeti_scale_entry_new_double(_("Vacancies:"),
                                                table, 3, &vvals.vacancies,
                                                0.0, 0.3, 0.001, 0.01, 3,
                                                G_CALLBACK(voronoi_g),
                                                drawable);

  /* visualization */
  vbox = gtk_vbox_new(FALSE, 0);
  gtk_container_set_border_width (GTK_CONTAINER (vbox), 4);
  gtk_notebook_append_page(GTK_NOTEBOOK(notebook), vbox,
                           gtk_label_new(_("Visualization")));

  frame = yeti_frame_new_in_box(_("Quantities Shown"), vbox, TRUE, TRUE);
  table = yeti_table_new_in_frame(3, 6, frame);

  for (i = 0; i < VORONOI_VALUE_LAST; i++) {
    vctrl.weights[i] = yeti_scale_entry_new_double(_(value_description[i]),
                                                   table, i, &vvals.weights[i],
                                                   0.0, 1.0, 0.01, 0.1, 2,
                                                   G_CALLBACK(voronoi_c),
                                                   drawable);
  }

  /* scale & shift */
  frame = yeti_frame_new_in_box(_("Scale & Shift"), vbox, TRUE, TRUE);
  table = yeti_table_new_in_frame(3, 2, frame);

  /* valscale */
  vctrl.valscale = yeti_scale_entry_new_double(_("R_escale (log):"),
                                               table, 4, &vvals.valscale,
                                               -1.5, 1.5, 0.01, 0.1, 2,
                                               G_CALLBACK(voronoi_c),
                                               drawable);

  /* valshift */
  vctrl.valshift = yeti_scale_entry_new_double(_("Shi_ft:"),
                                               table, 5, &vvals.valshift,
                                               -4.0, 4.0, 0.01, 0.1, 2,
                                               G_CALLBACK(voronoi_c),
                                               drawable);

  voronoi_celltype_update_dependants();
  voronoi_cellsize_changed(GTK_ADJUSTMENT(vctrl.cellsize));
  gtk_widget_show_all(dlg);
  in_param_reset = FALSE;
  voronoi(drawable, VORONOI_RECOMPUTE_GRID); /* preview image */

  ok = FALSE;
  do {
    run = gimp_dialog_run(GIMP_DIALOG(dlg));
    switch (run) {
      case VORONOI_RESPONSE_REVERT:
      voronoi_refresh_controls(&vvals_old);
      break;

      case VORONOI_RESPONSE_RESET:
      voronoi_refresh_controls(&vvals_defaults);
      break;

      case VORONOI_RESPONSE_ABOUT:
      about_dialog(dlg);
      break;

      default:
      ok = TRUE;
      break;
    }
  } while (!ok);

  gtk_widget_destroy(dlg);

  return run == GTK_RESPONSE_OK;
}

/*
 * The setup function.
 */

static void
voronoi(GimpDrawable *drawable,
        guint recompute)
{
  gint chmax;
  gdouble pcellsize;
  gint xwsq, xhsq, iter;
  gint step = 80;

  if (in_param_reset)
    return;

  if ((recompute & VORONOI_RECOMPUTE_GRID) && squares) {
    destroy_neighbourhoods(squares, wsq, hsq);
    destroy_grid_contents(squares, wsq, hsq);
    g_free(squares);
    squares = NULL;
  }

  init_voronoi(drawable);
  chmax = bpp;

  g_rand_set_seed(rng, vvals.seed);
  pcellsize = vvals.cellsize;
  if (preview_mode && vvals.cellsize > PREVIEW_CELL_SIZE_MAX)
      pcellsize = PREVIEW_CELL_SIZE_MAX;

  if ((recompute & VORONOI_RECOMPUTE_GRID)) {
    squares = random_grid(vvals.celltype, voronoi_width, voronoi_height,
                          pcellsize,
                          vvals.sigma, vvals.lambda,
                          vvals.htil, vvals.vtil,
                          vvals.interstit, vvals.vacancies,
                          &wsq, &hsq);
    xwsq = wsq + 2*SQBORDER;
    xhsq = hsq + 2*SQBORDER;
    progress += 0.1*max_progress;
    yeti_progress_update(progress, max_progress);
    for (iter = 0; iter < xwsq*xhsq; iter++) {
      find_voronoi_neighbours_iter(squares, wsq, hsq, iter);
      if (iter % step == 0) {
        progress += (0.4 * step)/(xwsq*xhsq)*max_progress;
        yeti_progress_update(progress, max_progress);
      }
    }
    /*
    if (vvals.celltype == VORONOI_GRID_RANDOM)
      remove_objects(squares, wsq, hsq, vvals.htil, vvals.vtil);
      */
    progress += (0.4 * (iter % step))/(xwsq*xhsq)*max_progress;
    yeti_progress_update(progress, max_progress);
  }
  commit(drawable, wsq, hsq, chmax);

  end_voronoi(drawable);
}

static void
init_voronoi(GimpDrawable *drawable)
{
  /* don't use time for seed when rendering preview in full size to get the
     same image as shown in preview */

  /* compute sizes and other stuff of this kind */
  alpha = -1000; /* to cause a segfault ;-> */
  has_alpha = gimp_drawable_has_alpha(drawable->drawable_id);
  if (preview_mode) {
    ix1 = iy1 = 0;
    ix2 = gdk_pixbuf_get_width(preview);
    iy2 = gdk_pixbuf_get_height(preview);
    bpp = gdk_pixbuf_get_n_channels(preview);
  }
  else {
    gimp_drawable_mask_bounds(drawable->drawable_id, &ix1, &iy1, &ix2, &iy2);
    bpp = gimp_drawable_bpp(drawable->drawable_id);
    if (has_alpha)
      alpha = bpp-1;
  }
  voronoi_width = ix2 - ix1;
  voronoi_height = iy2 - iy1;

  value_shift = vvals.valshift;
  value_scale = exp(vvals.valscale);

  max_progress = 2*(ix2 - ix1 + 1)*(iy2 - iy1 + 1);
  progress = 0;

  /* generate gradient, this is the right place to do it */
  if (vvals.cm == VORONOI_CM_GRADIENT) {
    if (preview_mode) {
      /* preview, done many times, don't regenerate what we already have */
      if (gradient == NULL) {
        gradient_p_width = CELL_SIZE_WIDTH;
        gradient = gimp_gradients_sample_uniform(gradient_p_width, FALSE);
      }
    }
    else {
      /* final rendering, done only once at the end */
      gimp_gradients_set_gradient(vvals.grad);
      g_free(gradient);
      gradient = gimp_gradients_sample_uniform(GRADIENT_SAMPLE_SIZE, FALSE);
    }
  }
}

static void
end_voronoi(GimpDrawable *drawable)
{
  if (preview_mode)
    gtk_widget_queue_draw(preview_image);
  else {
    gimp_drawable_flush(drawable);
    gimp_drawable_merge_shadow(drawable->drawable_id, TRUE);
    gimp_drawable_update(drawable->drawable_id,
                         ix1, iy1, (ix2 - ix1), (iy2 - iy1));
  }
}

/* give the abstract diagram some visual interpretation */
static void
commit(GimpDrawable *drawable,
       gint wsq,
       gint hsq,
       const gint chmax)
{
  gpointer pr;
  GimpPixelRgn dest_rgn;

  if (preview_mode) {
    dest_rgn.x = dest_rgn.y = 0;
    dest_rgn.w = dest_rgn.h = PREVIEW_SIZE;
    dest_rgn.bpp = gdk_pixbuf_get_n_channels(preview);
    dest_rgn.rowstride = gdk_pixbuf_get_rowstride(preview);
    dest_rgn.data = gdk_pixbuf_get_pixels(preview);
    commit_one_tile(&dest_rgn, wsq, hsq, chmax, gradient_p_width);
  }
  else {
    gimp_pixel_rgn_init(&dest_rgn, drawable,
                        ix1, iy1, (ix2 - ix1), (iy2 - iy1),
                        TRUE, TRUE);
    /* iterate by tiles */
    for (pr = gimp_pixel_rgns_register(1, &dest_rgn);
         pr != NULL;
         pr = gimp_pixel_rgns_process(pr)) {
      commit_one_tile(&dest_rgn, wsq, hsq, chmax, GRADIENT_SAMPLE_SIZE);
      progress += dest_rgn.w * dest_rgn.h;
      yeti_progress_update(progress, max_progress);
    }
  }
}

/* give the abstract diagram some visual interpretation (one tile) */
static void
commit_one_tile(GimpPixelRgn *rgn,
                gint wsq,
                gint hsq,
                const gint chmax,
                const gint gradwidth)
{
  static gint check_colors[] = {
    CHANNEL_MAX_VALUE*GIMP_CHECK_DARK,
    CHANNEL_MAX_VALUE*GIMP_CHECK_LIGHT,
  };

  VoronoiObject *owner, *line_start;
  VoronoiCoords z, zline, tmp;
  gint hsafe, vsafe;
  gint index = 0, channel, chc;
  gint x, y;
  guchar *p;
  gdouble r = 0, f, wsum;
  guint i;

  wsum = EPS;
  for (i = 0; i < VORONOI_VALUE_LAST; i++)
    wsum += vvals.weights[i];

  zline.x = wsq*(gdouble)rgn->x/voronoi_width;
  zline.y = hsq*(gdouble)rgn->y/voronoi_height;
  line_start = find_owner(squares, wsq, hsq, &zline);
  vsafe = 0;

  for (y = rgn->y; y < rgn->y + rgn->h; ) {
    hsafe = 0;
    z = zline;
    owner = line_start;

    neighbourize(owner->ne, &owner->pos);
    compute_segment_angles(owner->ne);

    tmp.y = zline.y;

    for (x = rgn->x; x < rgn->x + rgn->w; ) {
      /*fprintf(stderr, "x = %f, y = %f\n", z.x, z.y);*/
      /* draw the pixel */
      p = rgn->data + (y - rgn->y)*rgn->rowstride + (x - rgn->x)*rgn->bpp;
      r = 0.0;
      for (i = 0; i < VORONOI_VALUE_LAST; i++) {
        if (vvals.weights[i])
          r += vvals.weights[i]*value_function[i](&z, owner);
      }
      r = r/wsum*value_scale + value_shift;
      if (vvals.grad_invert)
        r = 1.0 - r;
      r = CLAMP(r, 0.0, 1.0);
      if (vvals.cm == VORONOI_CM_GRADIENT)
        index = r*(gradwidth - 1);
      for (channel = 0; channel < chmax; channel++) {
        if (vvals.cm == VORONOI_CM_GRADIENT) {
          p[channel] = CHANNEL_MAX_VALUE*gradient[4*index + channel];
          /* add checker in preview mode */
          if (preview_mode && has_alpha) {
            f = gradient[4*index + 3];
            chc = (x-ix1)%(2*GIMP_CHECK_SIZE)/GIMP_CHECK_SIZE
                  + (y-iy1)%(2*GIMP_CHECK_SIZE)/GIMP_CHECK_SIZE;
            p[channel] = p[channel]*f + (1-f)*check_colors[chc%2];
          }
        }
        else {
          if (channel == alpha)
            p[channel] = CHANNEL_MAX_VALUE;
          else
            p[channel] = CHANNEL_MAX_VALUE*r;
        }
      }

      /* move right */
      x++;
      if (hsafe-- == 0) {
        tmp.x = wsq*(gdouble)(x - ix1)/voronoi_width;
        owner = move_along_line(owner, &z, &tmp, &hsafe);
        neighbourize(owner->ne, &owner->pos);
        compute_segment_angles(owner->ne);
        z.x = tmp.x;
      }
      else
        z.x = wsq*(gdouble)(x - ix1)/voronoi_width;
    }

    /* move down */
    y++;
    if (vsafe-- == 0) {
      tmp.x = wsq*(gdouble)rgn->x/voronoi_width;
      tmp.y = hsq*(gdouble)(y - iy1)/voronoi_height;
      line_start = move_along_line(line_start, &zline, &tmp, &vsafe);
      zline.y = tmp.y;
    }
    else
      zline.y = hsq*(gdouble)(y - iy1)/voronoi_height;
  }
}

/* option menu callbacks */
static void
voronoi_celltype_update_dependants(void)
{
  gtk_widget_set_sensitive(vctrl.deform_frame,
                           vvals.celltype != VORONOI_GRID_RANDOM);
}

static void
voronoi_celltype_changed(GtkWidget *widget,
                         gpointer data)
{
  VoronoiGridType celltype;

  celltype = GPOINTER_TO_UINT(g_object_get_data(G_OBJECT(widget),
                                                "gimp-item-data"));
  if (vvals.celltype == celltype)
    return; /* no change, do nothing */
  vvals.celltype = celltype;
  voronoi_celltype_update_dependants();
  voronoi(drawable, VORONOI_RECOMPUTE_GRID);
}

static void
voronoi_cm_changed(GtkWidget *widget,
                   gpointer data)
{
  VoronoiColoringMethod cm;

  cm = GPOINTER_TO_UINT(g_object_get_data(G_OBJECT(widget), "gimp-item-data"));
  if (vvals.cm == cm)
    return; /* no change, do nothing */
  if (is_rgb) { /* this should be always true here, but... */
    gtk_widget_set_sensitive(vctrl.grad, cm == VORONOI_CM_GRADIENT);
    gtk_widget_set_sensitive(vctrl.grad_label, cm == VORONOI_CM_GRADIENT);
  }
  vvals.cm = cm;
  voronoi(drawable, VORONOI_RECOMPUTE_COLOR);
}

static void
voronoi_gradient_changed(const gchar *name,
                         gint width,
                         const gdouble *grad_data,
                         gboolean dialog_closing,
                         gpointer user_data)
{
  if (strncmp(vvals.grad, name, GRADIENT_NAME_SIZE-1) == 0)
    return;
  gimp_gradients_set_gradient(name);

  strncpy(vvals.grad, name, GRADIENT_NAME_SIZE-1);
  vvals.grad[GRADIENT_NAME_SIZE-1] = '\0';

  gradient = g_realloc(gradient, width*sizeof(gdouble));
  memcpy(gradient, grad_data, width*sizeof(gdouble));
  gradient_p_width = width/4;

  voronoi(drawable, VORONOI_RECOMPUTE_COLOR);
}

static void
voronoi_cellsize_changed(GtkAdjustment *adj)
{
  gdouble pcellsize;

  if (!preview_mode)
    return;

  pcellsize = vvals.cellsize;
  if (vvals.cellsize > PREVIEW_CELL_SIZE_MAX) {
    gchar *s;

    s = g_strdup_printf(_("Preview 1:%.1f"),
                        vvals.cellsize/PREVIEW_CELL_SIZE_MAX);
    gtk_frame_set_label(GTK_FRAME(vctrl.preview_frame), s);
    g_free(s);
  }
  else
    gtk_frame_set_label(GTK_FRAME(vctrl.preview_frame), _("Preview"));
}


/* though this is obscure, we use it to reset/revert all values _and_ refresh
   the dialog controls */
static void
voronoi_refresh_controls(const VoronoiValues *new_vvals)
{
  guint i;

  /* we cannot reset some field and don't want to reset some others, so copy
     them one by one  */
  in_param_reset = TRUE;
  vvals.htil = new_vvals->htil;
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vctrl.htil), vvals.htil);

  vvals.vtil = new_vvals->vtil;
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vctrl.vtil), vvals.vtil);

  vvals.celltype = new_vvals->celltype;
  gimp_option_menu_set_history(GTK_OPTION_MENU(vctrl.celltype),
                               GUINT_TO_POINTER(vvals.celltype));

  vvals.cellsize = new_vvals->cellsize;
  gtk_adjustment_set_value(GIMP_SCALE_ENTRY_SCALE_ADJ(vctrl.cellsize),
                           vvals.cellsize);

  vvals.sigma = new_vvals->sigma;
  gtk_adjustment_set_value(GIMP_SCALE_ENTRY_SCALE_ADJ(vctrl.sigma),
                           vvals.sigma);

  vvals.lambda = new_vvals->lambda;
  gtk_adjustment_set_value(GIMP_SCALE_ENTRY_SCALE_ADJ(vctrl.lambda),
                           vvals.lambda);

  voronoi_celltype_update_dependants();

  vvals.valscale = new_vvals->valscale;
  gtk_adjustment_set_value(GIMP_SCALE_ENTRY_SCALE_ADJ(vctrl.valscale),
                           vvals.valscale);

  vvals.valshift = new_vvals->valshift;
  gtk_adjustment_set_value(GIMP_SCALE_ENTRY_SCALE_ADJ(vctrl.valshift),
                           vvals.valshift);

  vvals.interstit = new_vvals->interstit;
  gtk_adjustment_set_value(GIMP_SCALE_ENTRY_SCALE_ADJ(vctrl.interstit),
                           vvals.interstit);

  vvals.vacancies = new_vvals->vacancies;
  gtk_adjustment_set_value(GIMP_SCALE_ENTRY_SCALE_ADJ(vctrl.vacancies),
                           vvals.vacancies);

  for (i = 0; i < VORONOI_VALUE_LAST; i++) {
    vvals.weights[i] = new_vvals->weights[i];
    gtk_adjustment_set_value(GIMP_SCALE_ENTRY_SCALE_ADJ(vctrl.weights[i]),
                             vvals.weights[i]);
  }

  vvals.cm = new_vvals->cm;
  if (!is_rgb)
    vvals.cm = VORONOI_CM_GRAY;
  gimp_option_menu_set_history(GTK_OPTION_MENU(vctrl.cm),
                               GUINT_TO_POINTER(vvals.cm));
  gtk_widget_set_sensitive(vctrl.grad, vvals.cm == VORONOI_CM_GRADIENT);
  gtk_widget_set_sensitive(vctrl.grad_label, vvals.cm == VORONOI_CM_GRADIENT);

  vvals.grad_invert = new_vvals->grad_invert;
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vctrl.grad_invert),
                               vvals.grad_invert);

  vvals.seed = new_vvals->seed;
  gtk_adjustment_set_value(GIMP_RANDOM_SEED_SPINBUTTON_ADJ(vctrl.seed),
                           vvals.seed);

  /* unable to update: grad -- FIXME */
  in_param_reset = FALSE;
  voronoi(drawable, VORONOI_RECOMPUTE_GRID);
}

/* GdkPixbuf RGB C-Source image dump {{{ */

#ifdef __SUNPRO_C
#pragma align 4 (voronoi_about_image)
#endif
#ifdef __GNUC__
static const guint8 voronoi_about_image[] __attribute__ ((__aligned__ (4))) = 
#else
static const guint8 voronoi_about_image[] = 
#endif
{ ""
  /* Pixbuf magic (0x47646b50) */
  "GdkP"
  /* length: header (24) + pixel_data (10800) */
  "\0\0*H"
  /* pixdata_type (0x1010001) */
  "\1\1\0\1"
  /* rowstride (180) */
  "\0\0\0\264"
  /* width (60) */
  "\0\0\0<"
  /* height (60) */
  "\0\0\0<"
  /* pixel_data: */
  "\342\224}\337\215x\325x]\330~e\224K!\1\0\0\0\0\0\0\0\0\0\0\0\225K!\337"
  "\215w\337\214w\227K\40\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\211L\"\321oR\374"
  "\325\227\325x]\235J\40G,\23\234J\40\316hH\360\265\212\343\226~\274G\36"
  "\276G\35\377\332J\377\364\313\377\376\376\377\376\376\377\364\311\377"
  "\332I\346\256\14\305\207\34\351\262\12\377\333O\377\365\317\377\376\376"
  "\377\376\376\377\354\241\377\323#\323\227\25\304\206\35\371\305\2\377"
  "\342n\377\374\362\377\376\376\377\367\332\377\335W\353\264\11\225M4|"
  "0@\327\235\23\377\334R\377\370\340\377\376\376\303R+\232J\40\205L\"\304"
  "R,\277I\37J-\24\0\0\0\0\0\0\0\0\0zK\"\322rU\357\265\212\275G\36\237J"
  "\40\212L!{K\"qF\37vH\40\304S-\357\264\212\335\210q\245I\37""0\35\15\0"
  "\0\0\24\14\5\211L\"\310[7\355\257\210\347\237\201\344\231\177\377\322"
  "\36\377\356\255\377\350\217\377\350\216\377\355\251\377\321\33\262p%"
  "e\25L\267v#\377\322\"\377\357\260\377\376\376\377\376\376\377\345\201"
  "\366\301\4\227P3s%E\321\225\26\377\332I\377\351\222\377\346\202\377\355"
  "\245\377\3251\301\202\36d\23L\246b+\377\315\6\377\344{\377\345\201\377"
  "\345\200^:\32\0\0\0\2\1\0\201M\"\313a\77\243J\40:$\20\0\0\0\0\0\0{L\""
  "\313cA\374\323\227\356\262\211\344\230~\332\203k\324uY\321nP\321nP\346"
  "\236\201\351\244\203\271H\36G,\23\0\0\0\0\0\0\0\0\0\0\0\0\203L\"\310"
  "[8\367\310\222\372\316\224\377\3266\377\333K\377\315\5\377\314\3\377"
  "\332F\377\3251\315\220\30\241\\.\321\225\26\377\3278\377\362\277\377"
  "\376\376\377\376\376\377\351\224\377\317\20\275~\40\253h)\351\262\12"
  "\377\337_\377\316\14\367\302\3\377\324(\377\333K\336\245\20\252g)\312"
  "\215\31\377\323'\377\330@\371\304\2\367\303\3nD\36\0\0\0\0\0\0\11\6\2"
  "\237J\40\325w[\243J\40S3\27H,\24\271H\36\350\244\203\347\240\201\322"
  "rU\323rU\326y^\334\206o\345\233\200\360\266\213\373\322\226\314dB|L\""
  "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\212L!\326y_\344\227~\311]:\377\344x\366"
  "\301\4\253h)\250e*\361\274\6\377\342q\377\323&\377\314\3\377\324)\377"
  "\346\204\377\375\367\377\376\376\377\376\376\377\366\323\377\340g\377"
  "\321\34\377\316\15\377\331D\377\320\24\271y\"\231R2\330\235\23\377\332"
  "I\377\330<\377\317\20\377\324*\377\343s\362\274\6\243^-\240[.\257I\37"
  "5!\16\0\0\0\0\0\0qF\37\310\\9\326y^\255I\37\246I\37\335\210r\360\265"
  "\212\305T/|M\"xJ!\177M\"\215L!\242J\40\276G\36\360\265\212\322rT\216"
  "K!W5\30\21\12\4\0\0\0\0\0\0\25\15\5\252I\37\347\237\201\304S-vI!\377"
  "\341j\347\260\13\213B9\207<;\342\252\16\377\337c\377\355\251\377\351"
  "\221\377\356\253\377\373\356\377\376\376\377\376\376\377\376\376\377"
  "\376\376\377\367\333\377\354\244\377\353\233\377\351\222\377\315\5\246"
  "a,v(C\311\214\32\377\330<\377\361\274\377\353\236\377\357\257\377\341"
  "k\351\262\12\220G7\214C8\330~d\236J\40)\31\13\0\0\0""8\"\17\257H\37\344"
  "\231\177\333\203k\326z_\375\326\230\322pR\212L!\0\0\0\0\0\0\0\0\0\0\0"
  "\0\0\0\0\227K\40\337\214w\344\227~\321oR\277H\37\235J\40\201M\"yJ!\204"
  "L\"\311^;\323rU\214L!\0\0\0\377\351\222\377\321\31\325\232\24\323\230"
  "\25\377\320\24\377\347\213\377\376\376\377\373\356\377\374\361\377\376"
  "\376\377\376\376\377\376\376\377\376\376\377\376\376\377\376\376\377"
  "\376\376\377\376\376\377\361\271\377\330=\346\256\14\321\225\26\374\311"
  "\1\377\341l\377\373\354\377\376\376\377\376\376\377\352\230\377\323#"
  "\336\245\20\335\243\20\326z`\324uX\237J\40A(\22""3\40\16\222K!\324uX"
  "\375\325\227\370\313\223\343\227~\253I\37$\26\12\0\0\0\0\0\0\0\0\0\0"
  "\0\0\0\0\0\206L\"\331\200f\335\210r\320mO\314cB\312_<\315gF\322pR\323"
  "rV\342\224}\260H\37)\31\13\0\0\0\377\371\342\377\346\204\377\332J\377"
  "\332H\377\345\200\377\370\335\377\351\222\377\337a\377\340e\377\353\236"
  "\377\357\261\377\351\223\377\355\247\377\371\345\377\376\376\377\376"
  "\376\377\376\376\377\376\376\377\354\241\377\336\\\377\332J\377\342r"
  "\377\363\306\377\376\376\377\376\376\377\376\376\377\373\354\377\350"
  "\220\377\335X\377\335W\261H\37\345\233\200\326z_\303P)\300I!\277H\37"
  "\305T.\357\264\212\367\307\222\312_<rF\37\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
  "\0\0\0\0\224K!\334\206o\330}d\252I\37sG\40^:\32c=\33\214L!\256I\37\342"
  "\224|\227K\40\0\0\0\0\0\0\377\363\307\377\336]\377\320\27\377\317\20"
  "\377\333N\377\351\222\377\322!\342\252\15\347\260\13\377\3250\377\325"
  "0\377\315\5\377\322\"\377\343w\377\372\350\377\376\376\377\376\376\377"
  "\365\316\377\337a\377\320\24\377\315\6\377\330>\377\354\242\377\376\376"
  "\377\376\376\377\376\376\377\360\266\377\332G\372\306\2\364\277\5\245"
  "I\37\330\177e\363\275\216\370\311\223\364\300\217\363\275\216\363\276"
  "\216\377\340\254\353\251\205\257I\37\33\21\7\0\0\0\0\0\0\0\0\0\0\0\0"
  "\0\0\0S3\27\302N'\354\253\206\330~e\330}c\252I\37O0\25\0\0\0\0\0\0\230"
  "K\40\342\224}\242J\40\3\1\0\17\11\4\377\346\205\377\314\2\267v#\257m"
  "'\362\274\6\377\337_\342\251\16\216E8\230R2\360\272\7\327\235\23\243"
  "^-\310\212\33\377\323'\377\356\255\377\376\372\377\375\367\377\351\221"
  "\377\316\13\271y\"\246b+\346\257\13\377\336]\377\371\345\377\364\313"
  "\377\370\335\377\345~\370\304\3\247c+\233U1\330~d\337\214w\306X3\310"
  "\\9\326y^\344\231\177\363\275\216\377\347\275\355\257\210\270H\36""2"
  "\36\15\0\0\0\0\0\0\0\0\0\0\0\0!\24\11\242J\40\336\213v\316hH\225K!\264"
  "H\36\316iJ\274G\36\213L!X6\30\232J\40\342\224}\266H\36\240J\40\255I\37"
  "\377\342r\354\265\11\215D8}1@\333\241\21\377\335Y\346\256\14\227P3\241"
  "[.\363\276\5\276\177\37g\27K\253g)\377\317\17\377\353\236\377\342n\377"
  "\341k\377\345\201\366\301\4\230Q2v)C\322\226\26\377\332J\377\337a\377"
  "\330=\377\335V\377\343w\361\273\6\227P3\207=;\315eD\242J\40iA\35W6\30"
  "\206L\"\246I\37\303Q*\364\277\216\377\336\245\323tW\206L\"\0\0\0\0\0"
  "\0\0\0\0\22\13\5\220K!\321nP\342\224}\245I\37\25\15\5""5\40\16~M\"\271"
  "H\36\325w\\\302N'\250I\37\345\232\177\340\216y\340\216z\313bA\377\330"
  ">\377\320\27\316\221\30\307\211\33\377\314\3\377\322\"\327\234\23\316"
  "\222\30\363\276\5\377\330=\356\267\10\304\206\34\341\250\16\377\330>"
  "\377\3265\355\266\10\351\262\12\377\324)\377\323'\327\234\23\311\214"
  "\32\375\311\0\377\336Z\360\272\7\310\213\32\345\255\14\377\331C\377\325"
  "-\342\252\16\333\242\21\201M\"\24\14\5\0\0\0\0\0\0\0\0\0\25\15\6\232"
  "J\40\337\215x\377\354\315\354\253\206\272H\36A(\22\0\0\0\24\14\5\212"
  "L!\312_<\361\270\214\313a\77nD\36\0\0\0\0\0\0\0\0\0Y7\30\260H\37\333"
  "\204l\341\222|\361\270\213\310\\8\243J\40pE\37\375\311\0\377\344|\377"
  "\330<\377\3266\377\341m\360\272\7\223K5}2@\331\236\22\377\335U\377\336"
  "\\\377\3277\377\334P\377\343t\360\272\7\233U1\223L5\345\255\14\377\337"
  "c\377\332H\377\330<\377\341l\377\324,\276\177\40j\33I\255j(\377\317\22"
  "\377\352\230\377\336\\\377\335V\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0g"
  "\77\34\313a\77\374\323\227\377\340\252\330}c\226K!\34\21\7\212L!\307"
  "Y5\354\253\206\347\237\201\256I\37$\26\12\0\0\0\0\0\0\0\0\0\0\0\0b<\33"
  "\263H\36\335\210r\340\217z\230K\40&\27\12\0\0\0\377\314\4\377\350\215"
  "\377\363\307\377\363\303\377\346\206\377\314\1\262q%\246b,\352\264\11"
  "\377\340e\377\370\337\377\363\305\377\366\326\377\343t\360\272\7\233"
  "U1\223K5\345\255\14\377\337c\377\366\322\377\364\311\377\363\307\377"
  "\330\77\324\231\25\237Y/\306\210\34\377\323&\377\356\255\377\372\350"
  "\377\371\343\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\"\25\11\263H\36\353"
  "\252\206\377\373\366\362\274\215\307Y4\216K!\307Y5\351\246\204\377\336"
  "\245\325w\\\214L!\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\240J\40\344\231"
  "\177\337\214v\222K!\0\0\0\0\0\0\377\333O\377\362\302\377\376\376\377"
  "\376\376\377\363\304\377\335W\377\316\14\377\314\3\377\331A\377\355\247"
  "\377\376\376\377\376\376\377\376\376\377\355\247\377\3264\354\266\11"
  "\351\261\12\377\324(\377\352\230\377\376\376\377\376\376\377\376\375"
  "\377\347\210\377\324*\376\313\0\377\321\35\377\343t\377\372\346\377\376"
  "\376\377\376\3763\37\16\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\223K!\334"
  "\207p\377\351\304\377\353\310\344\231\177\310\\9\351\245\204\377\347"
  "\276\363\274\215\305T/b<\33\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0+\32\13\266"
  "H\36\356\260\210\340\216y\226K!\0\0\0\0\0\0\377\354\240\377\376\376\377"
  "\376\376\377\376\376\377\373\355\377\350\216\377\334R\377\333M\377\345"
  "\201\377\367\333\377\376\376\377\376\376\377\376\376\377\375\371\377"
  "\353\233\377\337`\377\336Z\377\347\213\377\371\343\377\376\376\377\376"
  "\376\377\376\376\377\362\300\377\341l\377\331D\377\335W\377\353\236\377"
  "\376\376\377\376\376\377\376\376\271H\36\223K!tG\40Q2\26J.\24U4\27b<"
  "\33tG\40\316gG\376\330\231\377\376\376\377\332\233\374\324\227\377\347"
  "\275\377\356\322\345\232\177\306V1\240J\40qE\37(\30\13\0\0\0\0\0\0\11"
  "\6\2lB\36\307Z6\370\313\223\343\225}\237J\40sF\40U4\27\377\330>\377\361"
  "\272\377\376\376\377\376\376\377\353\235\377\323$\333\242\21\325\233"
  "\24\377\317\22\377\347\207\377\376\376\377\376\376\377\376\376\377\355"
  "\246\377\3250\345\255\14\336\245\20\377\321\31\377\347\213\377\376\376"
  "\377\376\376\377\372\350\377\340g\370\304\3\315\221\30\343\253\15\377"
  "\330<\377\361\271\377\376\376\377\376\376\332\202j\335\211r\317iJ\304"
  "R+\302O(\304R,\306W3\311^;\341\222|\377\334\240\377\376\376\377\376\376"
  "\377\376\376\377\376\376\377\376\376\377\347\276\365\303\220\342\223"
  "|\315gG\270H\36\230K\40}M\"\242J\40\311^;\344\231\177\377\341\254\363"
  "\274\215\341\220{\316iJ\311^:\377\315\10\377\352\226\377\347\213\377"
  "\347\211\377\343s\356\270\10\222J5\2049<\335\244\20\377\336[\377\355"
  "\247\377\347\213\377\354\240\377\344x\362\274\6\231R2\211\77:\336\244"
  "\20\377\336Z\377\352\230\377\347\211\377\356\254\377\3279\307\211\33"
  "r$E\246b,\377\315\7\377\351\224\377\343u\377\343u\221K!\301L$\334\206"
  "o\363\277\216\365\301\217\365\301\217\366\304\221\370\311\223\377\345"
  "\271\377\376\376\377\366\347\377\376\376\372\317\225\357\263\212\365"
  "\303\220\376\330\231\377\346\272\377\357\324\377\332\232\357\263\211"
  "\340\217{\324vZ\341\222|\370\312\223\377\353\312\365\301\217\363\276"
  "\216\373\321\226\377\334\237\375\326\230\377\321\31\377\332I\377\314"
  "\1\375\312\0\377\331A\371\305\2\253h)\241\\.\352\263\12\377\340f\377"
  "\323'\374\311\1\377\321\35\377\343u\370\303\3\245a,\230R2\345\255\14"
  "\377\337`\377\317\23\374\310\1\377\325/\377\332H\330\236\23\234V0\275"
  "}\40\377\320\30\377\3265\360\272\7\361\273\6\237J\40\323sV\323sW\305"
  "U/\327z`\351\246\204\374\323\227\377\355\316\377\362\334\375\326\230"
  "\366\305\221\377\376\376\364\277\217\305V1\302P(\314cB\327{a\344\230"
  "~\362\272\214\377\335\243\377\360\326\377\345\267\377\355\317\374\324"
  "\227\342\223|\306V1\277H\37\307Z6\353\251\205\332\202i\377\340h\366\301"
  "\4\251f*\245`,\357\271\7\377\333L\377\314\1\373\307\1\377\330<\377\327"
  ";\322\226\26\235W0\306\210\33\377\324(\377\331D\370\303\3\361\274\6\377"
  "\325/\377\320\26\274|!\236X0\333\241\21\377\333L\377\324+\372\307\2\377"
  "\317\22\377\340g\354\266\11\233U1\234V1\310[7\256I\37\223K!yJ!\222K!"
  "\272H\36\352\246\204\377\364\342\360\267\213\317jJ\341\220{\377\345\271"
  "\364\301\217\305U0`;\32d>\34\205L\"\237J\40\276G\36\336\212t\377\344"
  "\266\360\266\213\377\337\251\323sV\247I\37`;\32""4\40\16\232J\40\325"
  "x]\330~e\377\341m\351\262\12\216E7\207=;\341\251\16\377\337a\377\347"
  "\213\377\346\206\377\357\263\377\324*\274}\40k\33I\257l'\377\320\25\377"
  "\354\243\377\345~\377\344x\377\351\221\377\314\3\243_-r$E\311\213\32"
  "\377\330<\377\355\251\377\347\207\377\352\230\377\341i\350\261\13\222"
  "K5\223L5oD\37B(\22\26\16\6\0\0\0yK!\316gG\374\324\227\365\303\220\317"
  "kK\224K!\311]:\363\275\216\366\304\220\305V0]9\32\0\0\0\0\0\0E*\23\262"
  "H\36\344\231\177\355\257\210\334\206o\377\343\263\325x\\\205L\"\0\0\0"
  "\0\0\0`;\32\277H\37\330}d\377\352\226\377\322\36\331\236\22\326\233\24"
  "\377\320\25\377\347\213\377\376\376\377\376\376\377\367\330\377\335X"
  "\360\272\7\312\214\32\347\257\13\377\332F\377\363\305\377\376\376\377"
  "\376\376\377\360\265\377\3278\341\250\16\314\220\30\371\305\2\377\341"
  "i\377\372\352\377\376\376\377\376\376\377\352\232\377\324(\345\254\14"
  "\345\255\14\2\1\0\0\0\0\0\0\0\27\16\6\250I\37\344\231\177\375\326\230"
  "\325x\\\233J\40,\33\14\245I\37\336\212t\370\311\223\307Y4]9\31\0\0\0"
  "\0\0\0X6\30\275G\36\353\251\205\320mN\306X3\366\304\221\331\201h\213"
  "L!\0\0\0\0\0\0\15\10\3\224K!\331\177f\377\372\347\377\347\211\377\333"
  "O\377\333L\377\346\202\377\370\336\377\376\376\377\376\376\377\376\376"
  "\377\357\263\377\337b\377\330\77\377\335Y\377\355\245\377\376\376\377"
  "\376\376\377\376\376\377\376\376\377\353\233\377\334T\377\331B\377\341"
  "l\377\362\301\377\376\376\377\376\376\377\376\376\377\374\361\377\352"
  "\230\377\337b\377\337c\0\0\0\0\0\0\0\0\0{L\"\315fE\373\321\226\340\216"
  "y\251I\37""5\40\16\0\0\0qE\37\307Y4\364\301\217\311_<a<\33\0\0\0\0\0"
  "\0kB\35\305T.\325v[\246I\37\246I\37\344\231\177\340\216z\225K!\0\0\0"
  "\0\0\0\0\0\0Y7\30\307Y5\377\362\276\377\333N\376\313\0\366\301\4\377"
  "\3262\377\353\233\377\376\376\377\376\376\377\375\366\377\345\200\377"
  "\322\"\371\305\2\377\320\30\377\342r\377\372\346\377\376\376\377\376"
  "\376\377\357\263\377\331D\371\305\2\365\300\4\377\3266\377\354\242\377"
  "\376\376\377\376\376\377\376\376\377\360\265\377\331C\366\301\4\357\271"
  "\7S3\27\5\3\1>&\21\265H\36\350\243\203\355\257\210\300J!V5\27\0\0\0\0"
  "\0\0\31\17\7\243J\40\341\222|\317jJnD\36\0\0\0\0\0\0\177M\"\315fF\256"
  "I\37`;\32\200M\"\323sW\350\243\203\247I\37\1\1\0\0\0\0\0\0\0!\24\11\267"
  "H\36\377\346\206\376\312\0\254i(\233U1\343\253\15\377\336\\\377\372\347"
  "\377\367\332\377\362\301\377\3278\317\222\27\231S2\305\207\34\377\323"
  "'\377\357\257\377\364\310\377\365\321\377\344z\365\301\4\245a,\235W0"
  "\351\262\12\377\340f\377\367\327\377\363\303\377\367\333\377\345\177"
  "\370\304\3\245a,\226N4\302N&\232J\40\234J\40\331\200g\377\333\234\322"
  "qS\207L\"\0\0\0\0\0\0\0\0\0\0\0\0uH\40\316hH\327{a\204L\"\0\0\0\14\7"
  "\3\217K!\274G\36jA\35\0\0\0E*\23\304S,\363\276\216\277H\37""4\40\16\0"
  "\0\0\0\0\0*\32\13\257I\37\377\345}\365\300\4\232S2\2027>\330\236\23\377"
  "\334R\377\337a\377\333L\377\342r\377\324)\274}\40n\37G\261o&\377\320"
  "\27\377\341j\377\327:\377\332G\377\342r\355\267\10\224L4\211@:\340\250"
  "\16\377\336^\377\334Q\377\3265\377\335W\377\344{\364\277\5\233U1\211"
  "@:\360\267\213\334\205n\316iI\367\310\222\350\242\203\261H\37(\31\13"
  "\0\0\0\0\0\0\0\0\0\0\0\0-\34\14\273G\36\343\225}\234J\40\0\0\0f\77\34"
  "\275G\36\207L\"\12\6\2\0\0\0\25\15\6\262H\36\355\257\210\317jJtH\40^"
  ":\32\213L!\263H\36\323rU\377\323#\377\325/\341\251\16\327\235\23\377"
  "\317\17\377\331D\352\262\12\323\227\25\373\310\1\377\336Z\362\274\6\315"
  "\221\30\352\263\12\377\333K\373\310\1\310\213\33\327\235\23\377\324("
  "\377\324)\341\250\16\334\242\21\377\321\32\377\331C\343\252\15\303\204"
  "\35\351\262\12\377\334Q\377\3264\347\260\13\340\247\17\335\211s\377\333"
  "\235\365\302\220\377\333\236\321pR\202L\"\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
  "\0@'\21\241J\40\346\234\200\275G\36i@\35\254I\37\310[7lB\36\0\0\0\0\0"
  "\0\0\0\0\241J\40\345\233\200\341\220{\261H\36\305V0\326z_\332\202j\273"
  "G\36\365\301\4\377\345\200\377\336Z\377\334Q\377\345\201\377\316\13\252"
  "g*v(D\306\210\34\377\3277\377\340f\377\331D\377\337_\377\332I\321\225"
  "\26u'D\231R2\367\303\3\377\346\203\377\336[\377\335W\377\347\212\377"
  "\320\27\257m'd\23L\271y\"\377\323&\377\353\237\377\340d\377\336]\265"
  "H\36\337\215x\377\340\252\355\257\210\272H\36:$\20\0\0\0\0\0\0\0\0\0"
  "\20\11\4f\77\34\255I\37\326y_\335\211s\345\234\200\304S,\324tX\277H\37"
  "W6\30\0\0\0\0\0\0\0\0\0\227K\40\337\215y\362\274\215\347\240\202\357"
  "\264\212\317kL\236J\40D*\23\377\316\16\377\351\222\377\371\345\377\370"
  "\336\377\354\241\377\321\31\273{!\226O4\323\230\25\377\331B\377\364\314"
  "\377\366\323\377\371\344\377\336\\\346\256\14\246b+\272z!\377\316\16"
  "\377\351\223\377\372\346\377\371\343\377\360\265\377\325.\314\220\30"
  "\244_-\324\231\25\377\330<\377\363\303\377\373\357\377\372\352oD\36\317"
  "kL\377\337\251\362\272\214\272H\36&\27\12=&\21c=\33\207L\"\246I\37\310"
  "[8\322pR\272H\36\233J\40\311^;\355\256\207\346\234\200\262H\36D*\23\0"
  "\0\0\0\0\0\0\0\0\221K!\333\205m\377\343\264\353\251\205\312`=\223K!*"
  "\31\13\0\0\0\377\340e\377\365\321\377\376\376\377\376\376\377\367\327"
  "\377\340d\377\316\14\364\277\5\377\322\"\377\346\205\377\376\375\377"
  "\376\376\377\376\376\377\354\241\377\330>\377\315\7\377\320\26\377\337"
  "c\377\365\320\377\376\376\377\376\376\377\373\357\377\345~\377\323'\377"
  "\315\6\377\325.\377\347\211\377\376\374\377\376\376\377\376\376\237J"
  "\40\322pR\377\341\256\370\311\223\304R,\240J\40\264H\36\305U/\323sV\335"
  "\211r\274G\36\203L\"7\"\17\5\3\1\203L\"\307Y4\363\275\215\320mO\231J"
  "\40""6!\17\0\0\0*\32\13\263H\36\353\251\205\377\332\232\314cB\220K!'"
  "\30\12\26\15\6\12\6\2\377\352\226\377\376\374\377\376\376\377\376\376"
  "\377\376\376\377\354\243\377\336]\377\332I\377\342p\377\363\303\377\376"
  "\376\377\376\376\377\376\376\377\367\333\377\346\203\377\334S\377\336"
  "[\377\352\231\377\375\370\377\376\376\377\376\376\377\376\376\377\355"
  "\251\377\336\\\377\331A\377\337b\377\357\262\377\376\376\377\376\376"
  "\377\376\376\327\337\343\315\327\333\273\310\316\241\264\274\376\376"
  "\376\0""4J\0""4J\322\333\337\340\346\351\346\353\355\345\353\355\0""4"
  "J\0""4J\0""4J\0""4J\0""4J\0""4J\0""4J\0""4J/YkJo~^\177\215k\211\226q"
  "\216\232p\215\2315^o;ct:bs2\\m\"Obp\274Rp\274R\200\304f\315\347\302\300"
  "\341\263q\275T;~\40""1j\33""2l\34.b\31D\222%\203\305j\177\303dB\216$"
  "\35@\20\0\0\0\0\0\0!H\22D\221%|\302a\277\341\262\206\306mF\227'G\230"
  "'A\214$\\\263:\225\315~\303\343\267\260\332\237\262\333\242\5""8M\303"
  "\317\324\261\301\307\227\254\265\376\376\376\376\376\376\306\322\327"
  "\333\342\345\351\355\357\357\362\364\356\362\363\346\353\355\0""4J\0"
  "4J\0""4J\0""4J\0""4J\0""4J\0""4J\31H[4^oHn}Ux\206*Vh\77fvLp\200Ru\204"
  "Qt\203Hn}9ar:|\40I\234(\217\313x\337\360\331\307\344\273v\277Z=\203!"
  "\27""2\15\15\34\7#K\23H\233(\217\312w\216\312vH\233(#J\23\0\0\0\10\22"
  "\4-`\31R\256-\236\321\212\342\361\334\222\314{J\236(+]\30M\245*\215\312"
  "v\316\350\303\357\367\354\245\325\223b\266A\0""4J\0""4J\0""4J\77fw\0"
  "4J\0""4J\0""4J\335\344\347\352\356\360\361\364\365\360\363\364\350\355"
  "\356\331\340\344\376\376\376\0""4J\0""4J\0""4J\0""4J\243\266\276\274"
  "\311\317\315\327\334\0""4J\36L_:bsNr\201[}\213b\202\217`\201\216Xz\210"
  "In}e\267EO\251+\217\312w\336\357\330\273\337\254m\272N9z\37\24+\13\24"
  "+\13!G\22F\227'\212\310q\215\312vH\232'\"J\23\0\0\0\16\36\7""3n\34a\265"
  "A\261\332\241\334\356\325\215\311uH\232'8x\37h\270I\263\333\244\374\375"
  "\373\376\376\376\271\336\252i\271J\6""9N\0""4J\0""4J\0""4J\0""4J\0""4"
  "J\0""4J\0""4J\345\352\354\353\357\361\352\356\360\342\350\353\0""4J\376"
  "\376\376\376\376\376\376\376\376\376\376\376\376\376\376\233\257\270"
  "\263\303\311\305\321\326\4""7M'SeBiyWy\207d\204\221j\210\225i\210\224"
  "a\201\216Qu\203\263\333\243\225\315\177~\303c\312\346\277\236\321\212"
  "U\25714p\35:|\40""9{\40""3n\34\77\206\"v\277Y}\302bA\214$\35>\20\17\40"
  "\10\24*\13""4p\35c\266C\263\333\244\303\343\267y\300]@\211#\77\210#{"
  "\301`\312\346\277\376\376\376\376\376\376\273\337\255k\272L$Qd\27GZ\3"
  "7L\0""4J\0""4J\0""4J\0""4J\0""4J\0""4J\0""4J\0""4J\0""4J\0""4J\0""4J"
  "\0""4J\0""4J\376\376\3764^oAhxHm}\0""4J\5""8N(TfDjzXz\210e\205\222k\211"
  "\226j\211\225b\202\217Rv\204\376\376\376\332\356\323\250\326\226\246"
  "\325\223t\276WL\242*a\265@o\274Qo\273Q`\264>K\240)U\2571_\264=5r\35)"
  "X\26""4o\34""9{\37""9z\37W\2604\244\324\221\233\320\207W\26042l\34A\214"
  "$\177\303d\317\350\305\376\376\376\370\373\366\254\330\232_\264=<dt/"
  "Yk\33I]\3""7L\10:O\5""8N\0""4J\0""4J\0""4J\6""9N\14>S\13=R\2""6L\0""4"
  "J\0""4J\0""4JBhxVy\207c\203\220j\210\225i\207\224\0""4J\"Ob>evRv\204"
  "_\200\216e\205\222d\204\221\\}\213Mq\200\323\352\312\304\343\267\244"
  "\324\221w\277Za\265@\215\312v\256\331\235\277\341\262\276\341\261\254"
  "\330\233\212\310r]\263;D\222%:|\40L\242*a\265@n\273Pm\273O^\263<\206"
  "\306mh\270I\77\210#+]\30=\202!t\276W\302\342\265\376\376\376\321\351"
  "\310\215\311uK\241)Lp\200\77fv\23CW\37M`#Pc!Na\27FZ\6""9N\30G[%Rd+Vh"
  ")Ug!Na\21BV\0""4J\0""4J\\~\213q\216\232~\230\243\204\235\250\203\235"
  "\247{\226\241\25EY1[lEk{Rv\204Xz\211Wz\210Os\202@gw\203\305iw\277Z\\"
  "\263:`\264\77m\273O\223\315}\254\330\233\265\334\245\254\330\233\223"
  "\314|l\272NG\230'@\211#[\2629\200\304f\231\317\203\241\323\215\231\317"
  "\204\201\304g]\263;Z\2628<\201!L\243*l\272M\205\306l\243\324\220\323"
  "\352\312\236\322\212a\265A^\263<Ux\206Hm}(Tf3]n8ar5_p,Wi\33I]0Zl=duB"
  "iyAhx8ar)Tg\22BV\0""4Jp\215\231\204\235\250\222\250\261\230\255\266\227"
  "\254\265\217\246\257\177\231\244\35K^1[m>fvDjzCjy;ct+WiC\220%>\205\""
  "J\240)\213\311s{\301_M\245*^\263<e\267D]\263<M\245*=\204\":|\40a\265"
  "@H\232'D\222%N\246+Q\255-N\247+D\223&=\204\"^\263<X\2615\210\307o\264"
  "\334\244\301\342\264\260\332\237\264\334\244\300\341\263\230\317\203"
  "\244\324\221Wy\207Jo~5^oAhxEk{Ciy9ar(Tf@gwMq\200Sv\205Qu\204In}9ar\""
  "Ob\4""7M|\227\242\221\247\261\236\262\272\244\267\277\243\266\276\233"
  "\257\270\214\243\255\1""5K\26EY#Pc)Ug(Tf\37M`\20AU\35\77\20""2l\33]\263"
  ";\252\327\231y\300]>\206\"2k\33""5q\35""2k\33)Y\27(W\26K\241)\215\311"
  "uW\2614/f\32(W\26+]\30)W\26!F\22(W\26J\236)\203\305ia\265@p\274Rq\274"
  "Sc\266Bw\277[\242\324\217\300\341\263\315\350\303Qu\203(Tf;ctGm|Kp\177"
  "In~\77fw.YkIn~Vx\207\\}\213Z|\212Ru\204Biy+Vh\15\77S\202\233\246\226"
  "\254\264\243\266\276\252\273\302\251\272\302\240\264\274\221\247\261"
  "z\225\241\0""4J\0""4J\6""9N\5""8N\0""4J\247\271\301*[\27""7w\36j\271"
  "K\265\334\245i\271I8x\37\23*\12\17\40\10\14\33\7\16\37\10""3n\34_\264"
  ">\255\330\233d\266D4p\35\17\40\10\6\15\3\3\10\2\23(\12""7u\36d\267D\242"
  "\323\216U\2571:}\40:}\40""5r\35O\252,\227\316\202\336\357\327\307\345"
  "\273Ek{'Sf:bsFl{Jo\177Hm}>ev-XjKp\177Xz\210^\177\214\\~\213Tw\205Djz"
  "-Xj\17@U\0""4J\224\252\263\242\265\275\250\271\301\247\271\300\236\262"
  "\272\217\246\257x\224\237[|\212\0""4J\0""4J\232\257\267\254\275\304\266"
  "\305\313O\251+B\215$i\271I\221\314{N\250+,^\30\16\37\10\23*\12\23*\13"
  "\23)\12""8y\37m\272N\262\333\242c\266B4o\34\16\37\10\0\0\0\0\0\0\32""9"
  "\16\77\210#{\301_\261\333\241a\265A3n\34\25-\13""5q\35c\266C\262\333"
  "\242\376\376\376\344\362\3361[m\37M`2\\m>evBiy@gw6_p%QdFk{Sv\204Xz\211"
  "Wy\210Ns\202\77fv(Tf\12<Q:cs\213\243\255\230\255\266\237\262\272\236"
  "\261\272\225\253\264\206\237\251o\215\231Ru\204i\210\224\211\241\253"
  "\242\265\275\263\303\311\276\313\321\231\317\204y\300]Y\2616b\266B;\177"
  "\40'U\26""3m\34""9z\37""9z\37""3o\34""8x\37k\272L\237\322\213S\256/-"
  "b\31\24+\13\32""8\16\33;\17\35>\20C\217%\203\305h\260\332\237`\265\77"
  "3m\34\33;\17""8x\37l\272M\274\337\256\376\376\376\356\367\353\0""4J\20"
  "AU#Pb.Yk3]n0[l'Se\26EY9brFl|Lq\200Ko\177Biy2\\n\33J]Wy\207Vx\207Nr\201"
  ">fv\216\245\257\215\244\256\205\236\250v\222\235_\200\215Ahxj\210\225"
  "\211\241\253\242\265\275\264\303\312\276\314\321\243\324\220\206\306"
  "m]\263;P\253,8w\36J\236)^\264=m\273On\273P`\265\77K\241)\\\2629~\302"
  "cE\223&,^\30""8y\37@\210#A\213$<\201!@\210#{\301_\235\321\211Q\255,3"
  "m\34@\210#G\230'e\267E\264\334\245\376\376\376\344\362\337_\177\215l"
  "\212\227r\217\233\30G[\34K^\32I\\\20AU\0""4J&Re3\\n8ar7`q.Yk\37M`e\204"
  "\221k\211\226j\210\225b\202\217Rv\204<dtv\222\235n\213\230\350\355\357"
  "\335\344\347\313\325\332c\203\220\203\234\247\234\260\270\255\276\305"
  "\270\306\314U\2571H\232'7v\36U\2572\\\263:\211\310p\253\327\231\275\340"
  "\257\276\340\260\255\330\234\215\311u`\265\77Q\255,:}\40N\247+j\271J"
  "|\301`\177\303dr\275UX\2615e\267E|\302aD\221%X\2615y\300\\\213\311s\216"
  "\312w\234\320\207\344\362\336\307\345\274g\206\223t\220\234z\225\240"
  "y\224\240q\216\232a\201\217\232\257\267\220\246\260\13=Q\30G[\36L_\34"
  "K^\24DXd\204\221r\217\233x\224\237w\223\236o\214\230_\200\216In~+Vh\355"
  "\360\362\351\355\357\336\345\347\314\326\333\376\376\376u\221\235\216"
  "\245\257\240\263\273\252\273\302.c\31$N\24""6t\36Z\2627S\257/~\303d\236"
  "\321\211\256\331\235\255\331\234\234\320\207|\301`P\254,X\2615`\264>"
  "\220\313y\265\334\246\313\346\300\317\350\305\300\341\263\240\323\215"
  "t\276WP\254,k\272L\235\321\211\304\343\267\332\356\323\336\360\330\317"
  "\350\305\266\335\247\234\320\207g\206\223t\221\235{\226\241z\225\240"
  "q\216\232b\202\217\225\253\264\213\243\255z\225\240\0""4J\222\250\261"
  "\235\261\271Vx\207j\211\225x\223\237~\230\243}\227\242u\221\235e\205"
  "\222Os\2021[m\14>R\343\351\353\330\340\343\306\321\326\376\376\376\376"
  "\376\376y\224\240\213\242\254\225\253\264\11\24\5$N\24G\231'\205\306"
  "kb\265AD\223&Q\255,^\264=^\263<P\254,C\220%5q\35Y\2616V\2602\201\304"
  "f\240\322\214\260\332\237\256\331\236\234\321\210{\301`P\254,a\265@\210"
  "\307p\265\334\245\324\352\313\341\361\333\332\356\323\301\342\264\231"
  "\317\204\256\331\235a\201\216n\214\230t\220\234s\220\234k\211\226[}\213"
  "\26FZ\30H[\24DX\10:O\231\256\266\244\266\276Tw\206i\210\224v\222\236"
  "|\227\242{\226\241s\220\234d\204\221Mr\2010Zl\0""4J\0""4J\0""4J\0""4"
  "J\376\376\376\376\376\376\203\234\247\223\251\262\0""4J\12\26\5/e\32"
  "W\2603\244\324\221d\266D4p\35,_\30""2k\33""1j\33,^\30!H\22E\224&\201"
  "\304f]\263<E\225&S\256.`\264>_\264=P\254,C\220%G\231'\207\307nU\2571"
  "m\272N\206\306m\221\313z\213\311sw\277Z\200\303e\312\346\277Sv\205`\201"
  "\216f\206\223e\205\222]~\214*Uh3]n6_p1[l%Qd\21BV\15>S\34J^$Qdn\213\230"
  "t\220\234s\220\233k\211\226[}\213Ek{\36L_\37M`\32I\\\15\77S\0""4J\0""4"
  "J\202\234\246\231\256\267\251\273\302\0""4J\17\40\10""4p\35d\266C\245"
  "\325\222X\2615/f\32\12\27\6\14\32\6\14\32\6\7\17\4+\\\27P\253,\234\321"
  "\210j\271K7w\36-a\31""2l\34""2k\33,^\30,^\30Q\255,\237\322\213m\273O"
  ";~\40E\224&I\235(G\231'K\242)\224\315~\343\362\336\222\251\262Kp\177"
  "Ru\204Pt\203/Zk@gwIn~Lp\177Gl|;cs'Sf%Qd4^o<du>ev8`q*Vh\26FYLp\2003]n"
  "<du>ev8ar,Wi\30G[\0""4J\0""4J\251\272\302\271\307\315\302\316\323*[\27"
  "4o\34b\266B\207\307nH\233(%P\24\1\4\1\0\0\0\0\0\0\11\24\5/d\32X\2615"
  "\250\326\226i\271J7v\36\21&\11\24,\13\32""7\16\32""7\16/d\32X\2615\250"
  "\326\226x\300\\>\204\"\40E\21$M\23(U\26M\246+\231\317\204\352\365\345"
  "f\205\222]~\214\7:O&Re>evOs\202Xz\210Z|\212Ux\206In~6_p6_pEk{Mr\201O"
  "s\202In};ct'Se;csKo\177Sv\205Ux\206Ps\202Ciy/Yk\24DX\0""4J\261\301\310"
  "\301\316\323\312\325\331P\253,J\237)S\256/\\\263:7v\36\26""0\14\12\27"
  "\5\30""4\15!G\22%P\24-a\31T\2570\243\324\220Y\26171h\33.b\31""9{\37\77"
  "\207#\77\207#9{\37Q\256-\240\322\214t\276V<\200!.d\31/e\32*Z\27I\234"
  "(\216\312w\334\356\324y\224\240p\215\231\17@T.XjEk{Vy\207`\200\216b\202"
  "\217]~\214Qu\203>ev@gwOs\202Wz\210Xz\211Rv\204Ek{1[lKp\177[|\212d\203"
  "\220e\205\222`\200\216Sv\205\77fw$Qc\2""5Ku\221\235\210\240\252\224\251"
  "\263\236\322\212\215\311u_\264>N\246+6t\36\33;\17,^\30<\200!F\226'K\240"
  ")I\234(J\237)\204\305jH\231'>\204\"P\254,l\272M{\301_{\301_l\272MP\254"
  ",\210\307oa\265@M\244*W\2604Y\2616O\251+C\220%t\276W\274\337\256\205"
  "\236\250|\227\242l\212\226.YkFl{Wy\207`\200\216b\202\220]~\214Qu\204"
  ">fvBiyRu\204Z|\212[}\212Ux\206Hm}3]nTw\205d\203\221l\212\227n\214\230"
  "i\207\224\\}\213Hm}-Xj\13=Q\201\233\246\224\252\263\240\263\274b\266"
  "BL\242*9{\37Q\256-P\254,6s\35L\244*n\273O\210\307o\224\315}\217\313x"
  "{\301`Z\2628B\216$h\270H\226\316\200\270\336\251\313\346\300\313\346"
  "\300\270\336\252\226\316\200h\270Is\276V\224\315~\247\326\224\251\326"
  "\227\232\317\205|\302`R\256.\216\312w\212\242\254\201\233\245q\216\232"
  "'Sf\77fwPt\202Y{\211\\}\213Wy\207Ko\1777`q>evMq\201Ux\206Wy\207Pt\203"
  "Cjz/YkUx\206e\205\222n\214\230p\215\231j\211\225]~\214Jo~/Yk\0""4J\207"
  "\237\252\232\256\267\246\270\2775q\35)Y\27\34<\17""8x\37Q\255,Q\256-"
  "\210\307o\257\331\236\257\332\237\237\322\213\200\303eT\2571H\232'k\272"
  "L\242\324\217\325\353\314\355\366\351\361\370\356\340\360\332\276\340"
  "\260\216\312w\204\305jm\273O\237\322\213\306\344\272\335\357\326\341"
  "\361\333\301\342\265\220\313yp\274R\207\240\252~\231\244n\214\230\203"
  "\234\2471[mBhyKp\177Nr\201In}Qu\203Jo~2\\mAhxJo~Kp\177Ek{8`q#PcPt\203"
  "`\200\216i\207\224j\211\225e\204\221Xz\210#Pc\26EY\1""5K\0""4J\230\255"
  "\266\244\267\276\20#\11\6\15\3\5\13\3'U\26I\234(Z\2628Q\255,_\264=`\264"
  ">R\256.E\224&4o\34K\241)\216\312vi\271I\212\310r\235\321\211\241\323"
  "\215\224\315}w\277[v\277Y\271\336\252\216\312w[\2628{\301_\216\312v\221"
  "\313z\204\306kn\273P\205\306k~\230\243u\221\235e\204\221\210\240\252"
  "\222\250\262\225\253\264Lp\177Z{\211a\201\216`\201\216Y{\211Jo\177/Y"
  "k7`q8ar2\\m%Qd\204\235\250\216\245\257Sv\205\\}\2137`q\77fw@gw:bs,Wi"
  "\30G[\0""4J\217\246\260\233\257\270\0\0\0\0\0\0\21%\11""5r\35b\265Am"
  "\273O9z\37""2k\33""2l\34-`\31\"J\23""2k\33]\263;\253\327\231\204\305"
  "jH\232(P\252,Q\255,L\242*N\250+\230\317\202\342\361\334\253\330\232]"
  "\263;A\213$H\233(I\235(P\254,\231\317\203\250\326\226m\213\227d\204\221"
  "t\220\234\205\236\250\220\246\260>fvTw\205b\202\217i\207\224h\207\224"
  "a\201\217Rv\204=du\35K^\36L_\30G[\13=Q\205\236\250\217\246\257!Na7`q"
  "Gl|Os\202Os\202In~<dt'Se\13=R\0""4J\213\243\255\0\0\0\0\0\0\30""5\15"
  ">\204\"w\277[r\275T;~\40\25.\13\15\34\7\10\22\4\20$\11""6t\36h\270H\270"
  "\335\251\233\320\206O\251+*[\27+]\30""1i\33\\\263:\254\330\232\373\375"
  "\372\270\336\251h\270H6t\36#K\23""7v\36h\270I\266\335\247\275\340\257"
  "\252\273\302\255\275\304\250\272\301|\227\242\206\237\251\77gwUw\206"
  "c\203\220j\210\225i\210\225b\202\217Sv\205>ev!Na\0""4J\0""4Jn\214\230"
  "\177\231\244\211\241\253)Ug@gwOs\202Wy\207Xz\210Qu\204Djz/Zk\23CX\0""4"
  "J\206\237\251\0\0\0\0\0\0\33:\17@\212#~\302ch\270H7u\36\22&\12\0\0\0"
  "\0\0\0\17!\10""5q\35d\266D\263\333\243\243\324\220S\256/,_\30""8y\37"
  "E\225&_\264=\257\331\236\376\376\376\264\334\244d\267D5q\35!G\22;\177"
  "!s\275V\303\343\267\300\342\263\257\277\306\262\302\310\256\276\305k"
  "\211\226\0""4J9brNr\201\\~\214c\203\220c\203\220\\}\213Mq\2017`q\32I"
  "]\265\304\313\275\312\320a\201\216r\216\232\15>S*VhAgxPt\202Xz\210Y{"
  "\211Rv\204Ek{0Zl\24DX\241\264\274\231\255\266\0\0\0\2\5\1\30""4\15=\203"
  "\"v\277YP\254,-a\31\11\24\5\0\0\0\0\0\0\13\27\6-b\31R\256.\236\321\212"
  "\232\317\205N\250+E\225&d\266C\205\306k\227\316\202\241\323\215\350\364"
  "\343\236\322\212S\256.E\224&F\227'A\214$o\274Q\276\341\261\262\333\242"
  ":bs\260\300\307\0""4J\3""6L\20AU\26EY\24DXOs\202Vy\207Vx\207Nr\201@g"
  "w*Uh\15>S\0""4J\13=R\34K^'Se*Ug$Qc:csJo~Ru\204Rv\204Lq\200\77fv*Uh\16"
  "\77T\254\275\304\244\266\276\30""4\15&R\25/e\32""4p\35_\264>@\210#\37"
  "B\21\14\32\6\6\16\3\34<\17-a\31:}\40C\220%z\301^\201\304gF\227'u\276"
  "X\250\326\226\317\351\305\347\363\342\352\365\346\311\346\276\270\336"
  "\251\211\310p\210\307o\213\310r~\303cc\266B\250\326\226\223\314}Biyw"
  "\223\236\10;P\34K^)Ug/Zk.Yk&Re\26FYAhx:bsk\211\226f\205\222Z|\212\13"
  "=R#Pc5^o\77fwBiy>ev-Xj<duDjzEk{\77fv1[m\35K^\261\301\307\260\300\307"
  "\250\271\3019{\40J\236(X\2615a\265@G\231'M\246+\77\207#,^\30#L\23;\177"
  "\40O\252,m\273N\200\303e\204\305jx\300\\i\271J\253\330\232\347\363\342"
  "\362\371\357\335\357\326\270\335\251\230\317\203\317\350\305\302\342"
  "\265\327\354\317\333\356\323\313\346\300\253\327\231\203\305ih\270I"};
/* }}} */

static void
about_close(GtkWidget **dialog)
{
  gtk_widget_destroy(*dialog);
  *dialog = NULL;
}

static void
about_dialog(GtkWidget *parent)
{
  static GtkWidget *about = NULL;
  GtkWidget *vbox, *hbox, *widget;
  GdkPixbuf *pixbuf;

  if (about) {
    gtk_window_present(GTK_WINDOW(about));
    return;
  }

  about = gtk_dialog_new_with_buttons(_("About"), GTK_WINDOW(parent),
                                      GTK_DIALOG_NO_SEPARATOR
                                      | GTK_DIALOG_DESTROY_WITH_PARENT,
                                      GTK_STOCK_CLOSE, GTK_RESPONSE_CLOSE,
                                      NULL);
  gtk_container_set_border_width(GTK_CONTAINER(about), 6);
  gtk_window_set_resizable(GTK_WINDOW(about), FALSE);
  gtk_window_set_position(GTK_WINDOW(about), GTK_WIN_POS_CENTER_ON_PARENT);

  vbox = GTK_DIALOG(about)->vbox;
  gtk_box_set_spacing(GTK_BOX(vbox), 12);

  hbox = gtk_hbox_new(FALSE, 6);
  gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

  pixbuf = gdk_pixbuf_new_from_inline(-1, voronoi_about_image, FALSE, NULL);
  widget = gtk_image_new_from_pixbuf(pixbuf);
  g_object_unref(pixbuf);

  gtk_box_pack_start(GTK_BOX(hbox), widget, FALSE, FALSE, 0);
  gtk_misc_set_alignment(GTK_MISC(widget), 0.5, 0.0);

  widget = gtk_label_new(NULL);
  gtk_box_pack_start(GTK_BOX(hbox), widget, TRUE, TRUE, 0);
  gtk_label_set_markup
    (GTK_LABEL(widget),
     _("<big><b>Voronoi</b> " VERSION "</big>\n"
       "The ultimate GIMP pattern generator.\n"
       "\n"
       "By David Ne\304\215as (Yeti).\n"
       "E-mail: <i>yeti@physics.muni.cz</i>\n"
       "Web: <i>http://trific.ath.cx/software/gimp-plugins/voronoi/ </i>\n"));
  gtk_label_set_line_wrap(GTK_LABEL(widget), TRUE);
  gtk_label_set_selectable(GTK_LABEL(widget), TRUE);

  g_signal_connect_swapped(about, "delete_event",
                           G_CALLBACK(about_close), &about);
  g_signal_connect_swapped(about, "response",
                           G_CALLBACK(about_close), &about);
  gtk_widget_show_all(about);
}

