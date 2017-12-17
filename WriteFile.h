#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <stdarg.h>

#ifndef FALSE
#define FALSE 0
#endif /* !FALSE */

#ifndef TRUE
#define TRUE 1
#endif /* !TRUE */

/*----------------------------------------------------------------------------*/
/* fatal error, print a message to standard error and exit
*/
static void error(char * msg)
{
	fprintf(stderr, "error: %s\n", msg);
	exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*/
/* memory allocation, print an error and exit if fail
*/
static void * xmalloc(size_t size)
{
	void * p;
	if (size == 0) error("xmalloc input: zero size");
	p = malloc(size);
	if (p == NULL) error("out of memory");
	return p;
}

/* open file, print an error and exit if fail
*/
static FILE * xfopen(const char * path, const char * mode)
{
	FILE * f = fopen(path, mode);
	if (f == NULL)
	{
		fprintf(stderr, "error: unable to open file '%s'\n", path);
		exit(EXIT_FAILURE);
	}
	return f;
}

/* close file, print an error and exit if fail
*/
static int xfclose(FILE * f)
{
	if (fclose(f) == EOF) error("unable to close file");
	return 0;
}

/* skip white characters and comments in a PGM file
*/
static void skip_whites_and_comments(FILE * f)
{
	int c;
	do
	{
		while (isspace(c = getc(f))); /* skip spaces */
		if (c == '#') /* skip comments */
		while (c != '\n' && c != '\r' && c != EOF)
			c = getc(f);
	} while (c == '#' || isspace(c));
	if (c != EOF && ungetc(c, f) == EOF)
		error("unable to 'ungetc' while reading PGM file.");
}

/* read a number in ASCII from a PGM file
*/
static int get_num(FILE * f)
{
	int num, c;

	while (isspace(c = getc(f)));
	if (!isdigit(c)) error("corrupted PGM or PPM file.");
	num = c - '0';
	while (isdigit(c = getc(f))) num = 10 * num + c - '0';
	if (c != EOF && ungetc(c, f) == EOF)
		error("unable to 'ungetc' while reading PGM file.");

	return num;
}

/* read a PGM image file
*/
double * read_pgm_image(char * name, int * X, int * Y)
{
	FILE * f;
	int i, n, depth, bin = FALSE;
	double * image;

	/* open file */
	f = xfopen(name, "rb"); /* open to read as a binary file (b option). otherwise,
							in some systems, it may behave differently */

	/* read header */
	if (getc(f) != 'P') error("not a PGM file!");
	if ((n = getc(f)) == '2') bin = FALSE;
	else if (n == '5') bin = TRUE;
	else error("not a PGM file!");
	skip_whites_and_comments(f);
	*X = get_num(f);               /* X size */
	skip_whites_and_comments(f);
	*Y = get_num(f);               /* Y size */
	skip_whites_and_comments(f);
	depth = get_num(f);            /* pixel depth */
	if (depth < 0) error("pixel depth < 0, unrecognized PGM file");
	if (bin && depth > 255) error("pixel depth > 255, unrecognized PGM file");
	/* white before data */
	if (!isspace(getc(f))) error("corrupted PGM file.");

	/* get memory */
	image = (double *)xmalloc(*X * *Y * sizeof(double));

	/* read data */
	for (i = 0; i<(*X * *Y); i++)
		image[i] = (double)(bin ? getc(f) : get_num(f));

	/* close file */
	xfclose(f);

	/* return image */
	return image;
}

/*----------------------------------------------------------------------------*/
/* read a 2D ASC format file
*/
double * read_asc_file(char * name, int * X, int * Y)
{
	FILE * f;
	int i, n, Z, C;
	double val;
	double * image;

	/* open file */
	f = xfopen(name, "rb"); /* open to read as a binary file (b option). otherwise,
							in some systems, it may behave differently */

	/* read header */
	n = fscanf(f, "%d%*c%d%*c%d%*c%d", X, Y, &Z, &C);
	if (n != 4 || *X <= 0 || *Y <= 0 || Z <= 0 || C <= 0) error("invalid ASC file");

	/* only gray level images are handled */
	if (Z != 1 || C != 1) error("only single channel ASC files are handled");

	/* get memory */
	image = (double *)xmalloc(*X * *Y * Z * C * sizeof(double));

	/* read data */
	for (i = 0; i<(*X * *Y * Z * C); i++)
	{
		n = fscanf(f, "%lf", &val);
		if (n != 1) error("invalid ASC file");
		image[i] = val;
	}

	/* close file */
	xfclose(f);

	return image;
}

/*----------------------------------------------------------------------------*/
/* read an image from a file in ASC or PGM formats
*/
double * read_image(char * name, int * X, int * Y)
{
	int n = (int)strlen(name);
	char * ext = name + n - 4;

	if (n >= 4 && (strcmp(ext, ".asc") == 0 || strcmp(ext, ".ASC") == 0))
		return read_asc_file(name, X, Y);

	return read_pgm_image(name, X, Y);
}

/*----------------------------------------------------------------------------*/
/* write curves into a PDF file. the output is PDF version 1.4 as described in
"PDF Reference, third edition" by Adobe Systems Incorporated, 2001
*/
void write_curves_pdf(double * x, double * y, int * curve_limits, int M,
	char * filename, int X, int Y, double width)
{
	FILE * pdf;
	long start1, start2, start3, start4, start5, startxref, stream_len;
	int i, j, k;

	/* check input */
	if (filename == NULL) error("invalid filename in write_curves_pdf");
	if (M > 0 && (x == NULL || y == NULL || curve_limits == NULL))
		error("invalid curves data in write_curves_pdf");
	if (X <= 0 || Y <= 0) error("invalid image size in write_curves_pdf");

	/* open file */
	pdf = xfopen(filename, "wb"); /* open to write as a binary file (b option).
								  otherwise, in some systems,
								  it may behave differently */

	/* PDF header */
	fprintf(pdf, "%%PDF-1.4\n");
	/* The following PDF comment contains characters with ASCII codes greater
	than 128. This helps to classify the file as containing 8-bit binary data.
	See "PDF Reference" p.63. */
	fprintf(pdf, "%%%c%c%c%c\n", 0xe2, 0xe3, 0xcf, 0xd3);

	/* Catalog, Pages and Page objects */
	start1 = ftell(pdf);
	fprintf(pdf, "1 0 obj\n<</Type /Catalog /Pages 2 0 R>>\n");
	fprintf(pdf, "endobj\n");
	start2 = ftell(pdf);
	fprintf(pdf, "2 0 obj\n<</Type /Pages /Kids [3 0 R] /Count 1 ");
	fprintf(pdf, "/Resources <<>> /MediaBox [0 0 %d %d]>>\nendobj\n", X, Y);
	start3 = ftell(pdf);
	fprintf(pdf, "3 0 obj\n");
	fprintf(pdf, "<</Type /Page /Parent 2 0 R /Contents 4 0 R>>\n");
	fprintf(pdf, "endobj\n");

	/* Contents object - graphic contents */
	start4 = ftell(pdf);
	fprintf(pdf, "4 0 obj\n<</Length 5 0 R>>\n"); /* indirect length in obj 5 */
	fprintf(pdf, "stream\n");
	stream_len = ftell(pdf);
	fprintf(pdf, "%.4f w\n", width); /* set line width */
	for (k = 0; k<M; k++) /* write curves */
	{
		/* an offset of 0.5,0.5 is added to point coordinates so that the
		drawing has the right positioning when superposed on the image
		drawn to the same size. in that case, pixels are drawn as squares
		of size 1,1 and the coordinate of the detected edge points are
		relative to the center of those squares. thus the 0.5, 0.5 offset.
		*/

		/* initate chain */
		i = curve_limits[k];
		fprintf(pdf, "%.4f %.4f m\n", x[i] + 0.5, Y - y[i] - 0.5); /* first point */

		/* add remaining points of the curve */
		for (j = i + 1; j<curve_limits[k + 1]; j++)
			fprintf(pdf, "%.4f %.4f l\n", x[j] + 0.5, Y - y[j] - 0.5);

		/* if the curve is closed, market as such */
		j = curve_limits[k + 1] - 1;
		if (x[i] == x[j] && y[i] == y[j]) fprintf(pdf, "h\n");

		/* end curve - stroke! */
		fprintf(pdf, "S\n");
	}
	stream_len = ftell(pdf) - stream_len; /* store stream length */
	fprintf(pdf, "\r\nendstream\n"); /* EOL must be CRLF before endstream */
	fprintf(pdf, "endobj\n");

	/* Contents' stream length object - the use of this indirect object
	for the stream length allows to generate the PDF file in a single
	pass, specifying the stream¡¯s length only when its contents have
	been generated. See "PDF Reference" p.40. */
	start5 = ftell(pdf);
	fprintf(pdf, "5 0 obj\n%ld\nendobj\n", stream_len);

	/* PDF Cross-reference table */
	startxref = ftell(pdf);
	fprintf(pdf, "xref\n");
	fprintf(pdf, "0 6\n");
	fprintf(pdf, "0000000000 65535 f\r\n"); /* EOL must be CRLF in xref table */
	fprintf(pdf, "%010ld 00000 n\r\n", start1);
	fprintf(pdf, "%010ld 00000 n\r\n", start2);
	fprintf(pdf, "%010ld 00000 n\r\n", start3);
	fprintf(pdf, "%010ld 00000 n\r\n", start4);
	fprintf(pdf, "%010ld 00000 n\r\n", start5);

	/* PDF trailer */
	fprintf(pdf, "trailer <</Size 6 /Root 1 0 R>>\n");
	fprintf(pdf, "startxref\n");
	fprintf(pdf, "%ld\n", startxref);
	fprintf(pdf, "%%%%EOF\n");

	/* close file */
	xfclose(pdf);
}

/*----------------------------------------------------------------------------*/
/* write curves into a TXT file
*/
void write_curves_txt(double * x, double * y, int * curve_limits, int M,
	char * filename)
{
	FILE * txt;
	int i, k;

	/* check input */
	if (filename == NULL) error("invalid filename in write_curves_txt");
	if (M > 0 && (x == NULL || y == NULL || curve_limits == NULL))
		error("invalid curves data in write_curves_txt");

	/* open file */
	txt = xfopen(filename, "wb"); /* open to write as a binary file (b option).
								  otherwise, in some systems,
								  it may behave differently */

	/* write curves */
	for (k = 0; k<M; k++) /* write curves */
	{
		for (i = curve_limits[k]; i<curve_limits[k + 1]; i++)
			fprintf(txt, "%g %g\n", x[i], y[i]);
		fprintf(txt, "-1 -1\n"); /* end of chain */
	}

	/* close file */
	xfclose(txt);
}

/*----------------------------------------------------------------------------*/
/* write curves into a SVG file
*/
void write_curves_svg(double * x, double * y, int * curve_limits, int M,
	char * filename, int X, int Y, double width)
{
	FILE * svg;
	int i, k;

	/* check input */
	if (filename == NULL) error("invalid filename in write_curves_svg");
	if (M > 0 && (x == NULL || y == NULL || curve_limits == NULL))
		error("invalid curves data in write_curves_svg");
	if (X <= 0 || Y <= 0) error("invalid image size in write_curves_svg");

	/* open file */
	svg = xfopen(filename, "wb"); /* open to write as a binary file (b option).
								  otherwise, in some systems,
								  it may behave differently */

	/* write SVG header */
	fprintf(svg, "<?xml version=\"1.0\" standalone=\"no\"?>\n");
	fprintf(svg, "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n");
	fprintf(svg, " \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
	fprintf(svg, "<svg width=\"%dpx\" height=\"%dpx\" ", X, Y);
	fprintf(svg, "version=\"1.1\"\n xmlns=\"http://www.w3.org/2000/svg\" ");
	fprintf(svg, "xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n");

	/* write curves */
	for (k = 0; k<M; k++) /* write curves */
	{
		fprintf(svg, "<polyline stroke-width=\"%g\" ", width);
		fprintf(svg, "fill=\"none\" stroke=\"black\" points=\"");
		for (i = curve_limits[k]; i<curve_limits[k + 1]; i++)
			fprintf(svg, "%g,%g ", x[i], y[i]);
		fprintf(svg, "\"/>\n"); /* end of chain */
	}

	/* close SVG file */
	fprintf(svg, "</svg>\n");
	xfclose(svg);
}
