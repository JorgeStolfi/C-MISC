#define PROG_NAME "mlhtovoc"
#define PROG_DESC "convert the Melhoramentos vocabulary to \".voc\" format"

/* Last edited on 2024-12-25 09:51:30 by stolfi */

/* 
  Programa para o pre-processamento do vocabulario da
  Melhoramentos versão KLS-0.
  
  Criado por stolfi em 95-06-09, baseado no
  "trans.l" do tomasz de 95-06-08. 
  
  O vocabulario da Melhoramentos contém entradas no estilo abaixo:
  
          [<palavra>]#<categoria>[#<canonico>]
          
  onde <canonico> é a forma canônica de <palavra>.  Na versão KLS-0,
  apenas os verbos tem formas derivadas diferentes da forma canônica
  (que é o infinitivo impessoal).
  
  Entradas com <palavra> vazia significam que a forma verbal
  do <canonico> indicada pela <categoria> não existe.

  O programa gera, para cada entrada que tem <palavra> não-vazia,
  duas entradas ("analítica" e "sintética"):
  
          <palavra><a1><categoria><a2><nsufpal><sufcan>#

          <canonico><s1><categoria><s2><nsufcan><sufpal>#
          
  onde 

    <a1>, <a2>, <s1>, <s2> são delimitadores;
     
    <sufcan> e <sufpal> são as terminações distintivas
     de <canonico> e <paalvra>, respectivamente;
        
    <nsufpal> e <nsufcan> são os respectivos comprimentos.
   
  Mais precisamente, se <maxpref> é o maior prefixo comum entre 
  <palavra> e <canonico>, então <sufcanonico> e <sufpalavra> 
  são os respectivos sufixos, isto é,

      <canonico> = <maxpref><sufcanonico>
      <palavra>  = <maxpref><sufpalavra>
      
  Note que <sufcanonico> e/ou <sufpalavra> podem ser vazios.
  Os campos <nsufpal> <nsufcan> são dígitos hexadetricimais
  (0-9 a-z) que indicam o comprimento do outro sufixo:
   
      <nsufpal>  = |<sufpalavra>|
      <nsufcan>  = |<sufcanonico>|
      
  A versão corrente do programa supõe que estes comprimentos 
  são sempre menores ou iguais a 35.  
   
  Uso:
  
    mlhtovoc <entrada> <saida>
    
  onde <entrada> é o nome completo do arquivo de entrada,
  no formato melhoramentos, e <saida> é o nome *sem extensão*
  dos arquivos de saída.  O programa grava tres arquivos:
    
    <saida>.voa com as entradas analíticas, 
    <saida>.vos com as entradas sintéticas.
    <saida>.vms com as mensagens de erro.
     
AUTHOR
    J. Stolfi, DCC/Unicamp, june 1995.
*/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h> 

#define false 0
#define true 1 

typedef int bool_t;
typedef unsigned char byte;

/* Comprimentos maximos dos campos: */
#define MAXNPAL 35
#define MAXNCAT 35
#define MAXNLIN 255
  
/* Separadores de campos na saida --- entradas analiticas: */
#define SEPA1 ":"
#define SEPA2 ":"
  
/* Separadores de campos na saida --- entradas sinteticas: */
#define SEPS1 "::"
#define SEPS2 ":"

long num_linhas = 0;        /* Conta linhas, incluindo vazias */
long num_vazias = 0;        /* Conta linhas vazias */
long num_erros = 0;         /* Conta linhas com formato invalido */
long num_inexistentes = 0;  /* Conta formas flexionadas inexistentes */
long num_palavras = 0;      /* Conta palavras do vocabulário */
long num_analiticas = 0;    /* Conta entradas analíticas gravadas */
long num_sinteticas = 0;    /* Conta entradas sintéticas gravadas */

/* Max numero de erros: */
#define MAXERR 100000

#ifdef MSDOS

#include <process.h>

#define RMODE "rt"
#define WMODE "wt"

#else

#define RMODE "r"
#define WMODE "w"

#endif

#ifdef SUNOS

/* Prototype declarations that are missing in stdio.h: */

int printf(char *, ...);
int scanf(char *, ...);

int fread  (char *ptr, int size, int nitems, FILE *stream);
int fwrite (char *ptr, int size, int nitems, FILE *stream);

int fputc(char, FILE*);
int fputs(char *, FILE*);
int fprintf(FILE *, char *, ...);
int fscanf(FILE *, char *, ...);
int fclose(FILE *);
int fflush(FILE *);
int _flsbuf(unsigned char, FILE*);
int _filbuf(FILE *);

#endif

/* Internal prototypes: */

int main(int argc, char *argv[]);
FILE *openfile(char *name, char *mode);
bool_t le_linha(byte *lin, int maxlen, long *nlin);
void processa_linha(byte *lin, long nlin);
void mens_linha(char *msg, long num, byte *lin, long nlin);

void define_tabela(bool_t tb[], char *lets);
long cata(bool_t tb[], byte *buf, long nbuf);

char *txtcat(char *a, char *b);
void check (bool_t cond, char *msg);

/* Prototype declarations that are missing in stdio.h: */
/*
int printf(char *, ...);
int scanf(char *, ...);

int fread  (char *ptr, int size, int nitems, FILE *stream);
int fwrite (char *ptr, int size, int nitems, FILE *stream);

int fputc(char, FILE*);
int fputs(char *, FILE*);
int fprintf(FILE *, char *, ...);
int fscanf(FILE *, char *, ...);
int fclose(FILE *);
int fflush(FILE *);
int _flsbuf(unsigned char, FILE*);
int _filbuf(FILE *);
*/

static bool_t tb_pal[256]; /* Bytes que podem aparecer nas palavras */
static bool_t tb_cat[256]; /* Bytes que podem aparecer nas categorias */ 

static char digit[36] = "0123456789abcdefghijklmnopqrstuvwxyz";

static FILE *mlh, *voa, *vos, *vms;

int main(int argc, char *argv[])
  { byte *lin = (byte*) malloc(MAXNLIN + 1);
    long nlin;
    
    check(argc == 3, "usage: mlhtovoc <mlhin> <outname>");
    mlh = openfile(argv[1], RMODE);
    voa = openfile(txtcat(argv[2], ".voa"), WMODE);
    vos = openfile(txtcat(argv[2], ".vos"), WMODE);
    vms = openfile(txtcat(argv[2], ".vms"), WMODE); 
    if (vms == NULL) vms = stderr;
    
    define_tabela(tb_pal, "abcdefghijklmnopqrstuvwxyzáàãâéêíóôõúüç-");
    define_tabela(tb_cat, "abcdefghijklmnopqrstuvwxyz0123456789");
    
    while (le_linha(lin, MAXNLIN, &nlin))
      { num_linhas++;
        processa_linha(lin, nlin); 
      }
   
    fprintf(vms, "\n");
    fprintf(vms, "%8ld linhas no dicionario de entrada\n", num_linhas);
    fprintf(vms, "\n");
    fprintf(vms, "%8ld linhas vazias\n", num_vazias);
    fprintf(vms, "%8ld linhas de flexões inexistentes\n", num_inexistentes);
    fprintf(vms, "%8ld linhas malformadas\n", num_erros);
    fprintf(vms, "\n");
    fprintf(vms, "%8ld palavras válidas no vocabulário\n", num_palavras);
    fprintf(vms, "%8ld entradas analítcas gravadas\n", num_analiticas);
    fprintf(vms, "%8ld entradas sintéticas gravadas\n", num_sinteticas);
    return(0);
  } 
  
FILE *openfile(char *name, char *mode)
  { FILE *f; 
    if (vms != NULL) fprintf(vms, "opening %s  mode = %s\n", name, mode);
    if (vms != stderr) fprintf(stderr, "opening %s  mode = %s\n", name, mode);
    f = fopen(name, mode);
    check(f != NULL, txtcat("unable to open ", name));
    return(f);
  }
 

bool_t le_linha(byte *lin, int maxlen, long *nlin)
  /* 
    Lê proxima linha de "mlh", terminada por '\n'
    Devolve em "*nlin" o comprimento real da linha 
    (sem contar o '\n' final).  Coloca os bytes lidos
    em *lin (até o máximo de "maxlen" bytes, ignorando o excesso).
    Dá erro se a última linha não tiver '\n' no fim.
    Resultado é "false" sse "stdin" está esgotado.
  */
  { int c = getc(mlh);
    byte *p = lin;
    long n = 0;
    if (c == EOF)
      { check(ferror(mlh) == 0, "erro de leitura em stdin"); 
        *nlin = 0;
        return(false);
      }
    else
      { while((c != EOF) && (c != '\n'))
      { if (n < maxlen) *p = (byte) c; 
        n++; p++;
            c = getc(mlh);
      }
    if (c == EOF)
      { check(ferror(mlh) == 0, "erro de leitura em stdin"); 
        check(false, "ultima linha de stdin não termina com NEWLINE");
      }
    *nlin = n;
    return(true);
      }
  }
  
void mens_linha(char *msg, long num, byte *lin, long nlin)
  /*
    Imprime uma mensagem "msg" referente à linha numero "num",
    cujo conteudo (truncado em MAXNLIN bytes) é "*lin" e 
    cujo comprimento é "*nlin".
  */
  { fprintf(vms, "** linha %8ld", num);
    if (lin != NULL)
      { long n = nlin; 
        char *ell;
        fprintf(vms, " = \"");
        if (n > MAXNLIN) n = MAXNLIN;
        if (n > 40) { n = 40-3; ell = "..."; } else { ell = ""; }
        while (n > 0) 
          { if (((*lin) >= 32) && ((*lin) != 127))
              { fputc((char)(*lin), vms); }
            else
              { fputc((char)(*lin), vms); }
            n--; lin++;
          }
        fprintf(vms, "%s\"", ell);
      }
    fprintf(vms, ": %s\n", msg);
    check(num_erros < MAXERR, "excesso de erros");
  }

void processa_linha(byte *lin, long nlin)
  /* 
    Verifica, analisa, e processa uma linha,
    dados o comprimento total "nlin" e os primeiros MAXNLIN
    caracteres "lin[0..min(nlin,MAXNLIN)-1]".  Ambos excluem
    o NEWLINE final.
  */
  { byte *pal, *cat, *can;
    byte *sufpal, *sufcan;
    long npal, ncat, ncan, nsufpal, nsufcan;
    
    /* Verifica tamanho e caracteres: */
    if (nlin > MAXNLIN) 
      { mens_linha("muito comprida (ignorada)", num_linhas, lin, nlin);
        num_erros++;
        return;
      }
    if (nlin == 0) 
      { mens_linha("vazia (ignorada)", num_linhas, lin, nlin);
        num_vazias++;
        return;
      }

    /* Divide em campos (pal,npal) (cat,ncat) (can,ncan) */
    { byte *r = lin; long n = nlin;
    
      /* Palavra: */
      pal = r; 
      npal = cata(tb_pal, pal, n); 
      r += npal; n -= npal;
      if (npal > MAXNPAL) 
        { mens_linha("palavra muito longa (ignorada)", num_linhas, lin, nlin);
          num_erros++;
          return;
        }

      /* Separador: */
      if ((n == 0) || (((char)*r) != '#')) 
        { mens_linha("falta primeiro '#' (ignorada)", num_linhas, lin, nlin);
          num_erros++;
          return;
        }
      else
        { r++; n--; }

      /* Categoria */
      cat = r; 
      ncat = cata(tb_cat, cat, n); 
      r += ncat; n -= ncat;
      if (ncat > MAXNCAT) 
        { mens_linha("categoria muito longa (ignorada)", num_linhas, lin, nlin);
          num_erros++;
          return;
        }

      /* Separador: */
      if ((n == 0) || (((char)*r) != '#')) 
        { mens_linha("falta segundo '#' (ignorada)", num_linhas, lin, nlin);
          num_erros++;
          return;
        }
      else
        { r++; n--; }

      /* Forma canônica: */
      can = r; 
      ncan = cata(tb_pal, can, n); 
      r += ncan; n -= ncan;
      if (ncan > MAXNPAL)
        { mens_linha("forma canônica muito longa (ignorada)", num_linhas, lin, nlin);
          num_erros++;
          return;
        }

      /* fim de linha: */
      if (n != 0)
        { mens_linha("linha mal formada (ignorada)", num_linhas, lin, nlin);
          num_erros++;
          return;
        }
    }

    /* Campos vazios: */
    if (ncan == 0)
      { mens_linha("forma canônica vazia", num_linhas, lin, nlin);
        num_erros++;
        return;
      }
    if (ncat == 0) 
      { mens_linha("categoria vazia (ignorada)", num_linhas, lin, nlin);
        num_erros++;
        return;
      }

    /* Determina maior prefixo comum: */
    sufpal = pal; nsufpal = npal;
    sufcan = can; nsufcan = ncan;
    while((nsufpal > 0) &&  (nsufcan > 0) && ((*sufpal)==(*sufcan)))
      { sufpal++; nsufpal--; sufcan++; nsufcan--; }

    if (nsufpal > 35)
      { mens_linha("sufixo da palavra muito longo", num_linhas, lin, nlin);
        num_erros++;
        return;
      }
    if (nsufcan > 35)
      { mens_linha("sufixo da forma canônica muito longo", num_linhas, lin, nlin);
        num_erros++;
        return;
      }

    if (npal == 0) 
      { mens_linha("entrada de forma inexistente", num_linhas, lin, nlin);
        num_inexistentes++;
      }
    
    /* Grava resultado: */
    num_palavras++;
    pal[npal] = 0; cat[ncat] = 0; can[ncan] = 0;
    fprintf(voa, "%s%s%s%s%c%s#\n", pal, SEPA1, cat, SEPA2, digit[nsufpal], sufcan);
    num_analiticas++;
    fprintf(vos, "%s%s%s%s%c%s#\n", can, SEPS1, cat, SEPS2, digit[nsufcan], sufpal);
    num_sinteticas++;
  }

void define_tabela(bool_t tb[], char *lets)
  /* Define "tb[c]=true" para todo "c" em "*lets"; demais "false". */
  { int k;
    for (k=0; k<256; k++) tb[k] = false;
    while ((*lets)!='\0') {tb[(byte)(*lets)] = true; lets++; }
  }

long cata(bool_t tb[], byte *buf, long nbuf)
  /* 
    Devolve o maior "k" entre 0 e "nbuf" tal que "buf[0..k-1]"; 
    está inteiramente no conjunto "tb".
  */
  { int k = 0;
    while ((nbuf > 0) && (tb[*buf])) { buf++; nbuf--; k++; }
    return(k); } 
    
char *txtcat(char *a, char *b)
  { char *r = malloc(strlen(a) + strlen(b) + 1);
    check(r != NULL, "memory exhausted");
    strcpy(r, a); strcat(r, b); 
    return(r);
  }

void check(bool_t cond, char *msg)
  { if(!cond) 
      { 
        if (vms) fprintf(vms, "** %s\n", msg); 
        if (vms != stderr) fprintf(stderr, "** %s\n", msg); 
        exit(1); 
      } 
  }
  
