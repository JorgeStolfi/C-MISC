#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

/* Undoes CGI-BIN formfield encoding
 *
 * Last edited 96-05-15 by stolfi
 *
 * NAME
 *    decode-url-form-data - decodes HTML form data
 *
 * SYNOPSIS
 *   echo "$QUERY_STRING" | decode-url-form-data FS PREFIX
 *
 * DESCRIPTION
 * 
 *   When an HTML form is sent, its field names and values
 *   are encoded as follows:
 *   
 *     [A-Za-z0-9] = unchanged
 *     [@._*-;]    = either unchanged or recoded as %HH
 *     SP          = +                        
 *     other       = recoded as %HH 
 *  
 *   Then the name and value of each field are listed, separated
 *   by "=". 
 *
 *   The name-value pairs of all fields are then concatenated, using
 *   either ";", "|", or "&" as the field separator. The HTML 2.0
 *   standard recommends ";", and allows "&"; but I have seen "|"
 *   used, too.  Does it depends on the form's METHOD?  Anyway, this
 *   program let's you specify the field separator character, as the
 *   first argument FS.
 *  
 *   This program inverts the above mapping, and writes each field
 *   value to a file whose name is the PREFIX concatenated
 *   with the field's name.  The PREFIX may be empty, or 
 *   contain "/"s (in which case the named directories must exist).
 *   If any of those files already exists, the program aborts.
 *  
 *   Note that an encoded field name or value must be a string 
 *   in the language
 *
 *     ( [a-zA-Z0-9@._*;] + { "+" } + [%][0-9A-F][0-9A-F] )*
 *
 *   It is a fatal error if the input does not fit into this format.
 *  
 *   We don't do anything smart about newlines; "%0D" becomes CR, 
 *   "%0A" becomes LF, and that is all. (The HTML 2.0 standard
 *   says they are encoded as CR-LF ("%0D%0A") pairs. That rule may be
 *   server-friendly, but it is not invertible. And it deosn't seem
 *   to be followed by some browsers... Argh!...)
 *
 *   Note that the files created by this program do not end with 
 *   newline, unless the field value already ended with "%0A".
 *
 *   The field names must be at most than MAXNAMELEN bytes long, begin
 *   with alphameric, and contain only alphamerics or {'.', '-', '_'}.
 *   There must be at most MAXFIELDS fields.  Any violation of these
 *   rules is a fatal error.
 * 
 * AUTHOR
 *   J. Stolfi, IC-UNICAMP, May 1996
 * */
 
#define MAXNAMELEN 255
#define MAXFIELDS 1024

/* PROTOTYPES */

int main(int argc, char **argv);
  /* 
    Main program */

char *get_prefix(int argc, char **argv);
  /* 
    Parses the command line options, returns the filename prefix */

char get_field_separator(int argc, char **argv);
  /* 
    Parses the command line options, returns the field
    separator character */

char *get_field_name(void);
  /* 
    Parses a field name from stdin.  Checks if the name is
    well-formed, and at most MAXNAMELEN bytes long.  Returns the
    decoded name if sucessful, else aborts with an error message.
  */

FILE *open_file(char *prefix, char *fname);
  /* 
    Opens the file for field "fname", checks for failure */ 
    
void write_field_value(FILE *f);
  /*
    Reads the field value from stdin, decodes it,
    and writes the result to file "f".
    Aborts in case of syntax error. */

char read_hex_code(void);
  /*
    Reads two hex digits and returns the equivalent
    character. */

int is_delim(char c);
  /* 
    Tests if "c" is a field/value delimiter:  {FS, '=', EOF} */
  
int is_hex_digit(char c);
  /*
    Tests if "c" is a uppercase hex digit [0-9A-Z]. */

int is_alphanum(char c);
  /*
    Tests if "c" is an ASCII letter or digit. */

int is_plain(char c);
  /*
    Tests if "c" is a field value character that may be 
    left unencoded. */

char *txtcat (const char *a, const char *b);
  /* 
    Concatenation of "a" and "b", allocated. */

void error(char *msg);
  /*
    Writes error message "msg" to "stderr" and aborts. */

/* IMPLEMENTATIONS */

#define unhex(x) \
  ( ((x) >= '0') && ((x) <= '9') ? (x) - '0' : \
  ( ((x) >= 'A') && ((x) <= 'F') ? (x) - 'A' + 10 : \
    -1 ))
    
static long num_read = 0;

static char FS; /* Field separator */

int main(int argc, char **argv)
  { long num_fields = 0;
    char *prefix;
    char c;
    /* Command line parsing: */
    if (argc != 3) 
      { error("bad arguments: FS PREFIX expected"); }
    FS = get_field_separator(argc, argv);
    prefix = get_prefix(argc, argv);
    /* Field loop: */
    while((c = getchar()) != EOF)
      { char *fname;
        FILE *vf; 
        num_read++;
        if (num_fields >= MAXFIELDS)
          { error("too many fields"); }
        ungetc(c, stdin); num_read--;
        fname = get_field_name();
        num_fields++;
        vf = open_file(prefix, fname);
        c = getchar(); num_read++;
        if (c == '=')
          { /* Writes value into file */ 
            write_field_value(vf);
            c = getchar(); num_read++;
          }
        free(fname);
        fclose(vf);
        if ((c != FS) && (c != EOF))
          { error("missing field separator"); }
      }
    return(0);
  }

char get_field_separator(int argc, char **argv)
  { if (
      (argv[1][0] == '\0') || 
      (argv[1][1] != '\0') || 
      (argv[1][0] == '=') || 
      (argv[1][0] == '%') || 
      (argv[1][0] == '+') || 
      is_alphanum(argv[1][0])
    ) 
      { error("bad field separator"); }
    return(argv[1][0]);
  }

char *get_prefix(int argc, char **argv)
  { return(argv[2]); }

char *get_field_name(void)
  { char c;
    int num_chars = 0;
    char *fname = (char*)malloc(MAXNAMELEN+1);
    if (fname == NULL) error("out of memory");
    while (! is_delim(c = getchar()))
      { num_read++;
        if (num_chars >= MAXNAMELEN)
          { error("field name too long"); }
        if (c == '%') { c = read_hex_code(); }
        if (is_alphanum(c))
          { /* OK */ }
        else if (( c == '.' ) || ( c == '_' ) || ( c == '-' ))
          { if (num_chars == 0)
              { error("field name must begin with alphanum"); }
          }
        else
          { error("bad character in field name"); }
        fname[num_chars] = c;
        num_chars++;
      }
    if (c != EOF) ungetc(c, stdin);
    fname[num_chars] = '\0';
    return(fname);
  }
      
int is_delim(char c)
  {
    return((c == FS) || (c == '=') || (c == EOF));
  }
  
int is_hex_digit(char c)
  {
    return(
      ((c >= '0') && (c <= '9')) || 
      ((c >= 'A') && (c <= 'F'))
    ); 
  }
  
char read_hex_code(void)
  { char x, y;
    int xv, yv;
    x = getchar(); num_read++;
    if (x == EOF)
      { error("unexpected EOF in hex code"); }
    xv = unhex(x);
    if (xv < 0) 
      { error("bad hex digit"); }
    y = getchar(); num_read++;
    if (y == EOF)
      { error("unexpected EOF in hex code"); }
    yv = unhex(y);
    if (yv < 0) 
      { error("bad hex digit"); }
    return((char)((xv<<8)+(yv)));
  }

FILE *open_file(char *prefix, char *fname)
  { char *path = txtcat(prefix, fname);
    FILE *f;
    struct stat s;
    if (lstat(path, &s) == 0) { error(txtcat("file already exists: ", path)); }
    f = fopen(path, "w");
    if (f == NULL) { error(txtcat("cannot create file ", path)); }
    free(path);
    return(f);
  }
    
void write_field_value(FILE *f)
  { char c;
    while (! is_delim(c = getchar()))
      { num_read++;
        if (c == '+')
          { fputc(' ', f); }
        else if (c == '%')
          { fputc(read_hex_code(), f); }
        else if (! is_plain(c))
          { error("bad character in field value"); }
        else
          { fputc(c, f); }
      }
    if (c != EOF) ungetc(c, stdin);
  }
    
int is_alphanum(char c)
  { return(
      ( (c >= 'a') && (c <= 'z') ) |
      ( (c >= 'A') && (c <= 'Z') ) |
      ( (c >= '0') && (c <= '9') )
    );
  }
        
int is_plain(char c)
  { return(
      is_alphanum(c) |
      (c != '_') |
      (c != '@') |
      (c != '*') |
      (c != '.') |
      (c != '-')
    );
  }
        
void error(char *msg)
  {
    fprintf(stderr, "decode-url-form-data: byte %ld: %s\n", num_read, msg);
    exit(1);
  }

char *txtcat (const char *a, const char *b)
  {
    char *r = malloc(strlen(a)+strlen(b)+1);
    if (r == NULL) { error("memory exhausted"); }
    (void) strcpy(r, a);
    (void) strcat(r, b);
    return(r);
  }

