/* 
  Last edited on 2009-02-10 08:30:43 by stolfi
  
  Programa para pre-processamento do dicion�rio de flex�es de S. Carlos.
  
  UNDER DEVELOPMENT
  
  Criado por stolfi em 96-01-10, baseado no
  "trans.l" do tomasz de 95-06-08. 
  
  USO
  
    dfltovfl <entrada> <saida>
    
  onde <entrada> � o nome completo do dicion�rio de flex�es de
  S. Carlos, no formato de qu�druplas, e <saida> � o nome *sem extens�o*
  dos arquivos de sa�da.  O programa grava dois arquivos:
    
    <saida>.vfa com o dicion�rio de flex�es anal�tico, 
    <saida>.vfs com o dicion�rio de flex�es sint�tico,
    <saida>.vms com as mensagens de erro.

  FUN��O
 
  O dicion�rio de S. Carlos, depois de suficientemente 
  massageado, produz um conjunto de quu�ntuplas ("entradas") 
  da forma
  
    (<palavra>,<ordem>,<classe>,<flexao>,<canonica>)
          
  onde  
    
    <canonica> � um string de at� MAXLPAL caracteres do alfabeto,
        que � uma palavra da l�ngua em sua flex�o can�nica.
  
    <classe> � uma letra que indica uma fun��o gramatical
        geral que a <canonica> pode ter, como verbo, 
        substantivo, etc. (Vide tabela adiante)
      
    <flexao> � um string de at� MAXNCAT-1 letras e d�gitos que
        indicam uma flex�o aplic�vel � <canonica> quando ela
        occorre como <classe>: g�nero e n�mero para substantivos,
        tempo e pessoa para verbos, etc.
       
    <palavra> � um string de at� MAXLPAL caracteres do alfabeto,
        que representa a forma da <canonica> flexionada conforme
        indicado por <classe> e <flexao>. 
        
    <ordem> � um d�gito de ordem que indica a import�ncia
        relativa das v�rias entradas de uma mesma <palavra>.
 
  A concatena��o de <classe> e <flexao> � chamada por n�s de
  <categoria>.
       
  A <palavra> pode ser vazia, significando que essa particular
  <flexao> da forma <canonica> definidamente n�o existe. 
  
  O "vocabul�rio" � o conjunto de todas as <palavras> que aparecem
  neste arquivo.
  
 O dicion�rio de flex�es ser� usado basicamente de tr�s 
  maneiras: 
  
    busca "existencial": dada uma m�scara,
    enumerar as <palavras> do vocabul�rio que satisfazem
    essa m�scara. (Em particular, dada uma �nica <palavra>,
    verificar se ela est� no vocabul�rio,)
  
    busca "sint�tica": dada uma tripla <canonica>, <classe>,
    <flexao>, determinar a <palavra> correspondente.
    
    busca "anal�tica": dada uma <palavra> do vocabul�rio, 
    determinar todas as qu�druplas <ordem>,<canonica>,<classe>,<flexao>
    associadas a ela.  
    
  Para implementar essas buscas de maneira eficiente
  usando automatos, codificaremos o dicion�rio de flex�es 
  por dois conjuntos de cadeias, representando aproximadamente a
  fun��o flexionadora e sua inversa. Especificamente, para cada
  entrada que tem <palavra> n�o-vazia,  criamos
  duas cadeias ("sint�tica" e "anal�tica"):
  
    <palavra><a1><ordem><classe><flexao><a2><nsufpal><sufcanonica>#

    <canonica><s1><classe><flexao><s2><nsufcan><sufpalavra># 
          
  onde 

    <a1>, <a2>, <s1>, <s2> s�o delimitadores, 
        sendo que <s1> � distinto de <a1>;
     
    <sufcanonica> e <sufpalavra> s�o as termina��es distintivas
        de <canonica> e <palavra>, respectivamente;
        
    <nsufpal> e <nsufcan> s�o os respectivos comprimentos.
   
  Mais precisamente, se <maxpref> � o maior prefixo comum entre 
  <palavra> e <canonica>, ent�o <sufcanonica> e <sufpalavra> 
  s�o os respectivos sufixos, isto �,

    <canonica> = <maxpref><sufcanonica>
    <palavra>  = <maxpref><sufpalavra>
      
  Note que <sufcanonica> e/ou <sufpalavra> podem ser vazios.
  Os campos <nsufpal> <nsufcan> s�o d�gitos hexatrisdecimais 
  (0-9a-z) que indicam o comprimento do outro sufixo:
   
    <nsufpal>  = |<sufpalavra>|
    <nsufcan>  = |<sufcanonica>|
      
  A vers�o corrente do programa sup�e que estes comprimentos 
  s�o sempre menores ou iguais a 35.  
   
*/

#define _GNU_SOURCE
#include <stdio.h> 
#include <string.h>
#include <stdlib.h>

#define false 0
#define true 1 

typedef int bool_t;
typedef unsigned char byte;

/* Comprimentos maximos dos campos: */
#define MAXLPAL 34
#define MAXLORD  1
#define MAXLCLS 34
#define MAXLFLX 34
#define MAXLLIN 255
  
/* Max numero de erros: */
#define MAXERR 10000 
  
/* Max numero de registros a processar: */
#define MAXREG 1000000
        
/* Separadores de campos -- dicion�rio sint�tico de flexoes: */
#define SEPS1 ":"
#define SEPS2 "=" 

/* Separadores de campos -- dicionario analitico de flexoes: */
#define SEPA1 "="
#define SEPA2 ":"
  
long num_linhas = 0;        /* Conta linhas, incluindo vazias */
long num_vazias = 0;        /* Conta linhas vazias */
long num_erros = 0;         /* Conta linhas com formato invalido */
long num_inexistentes = 0;  /* Conta formas flexionadas inexistentes */
long num_palavras = 0;      /* Conta palavras do vocabul�rio */
long num_analiticas = 0;    /* Conta entradas do dic anal�tico de flexoes */
long num_sinteticas = 0;    /* Conta entradas do dic sint�tico de flexoes */

/* Arquivos */
static FILE 
  *dfl,   /* Dicion�rio de S. Carlos, formato quintuplas */
  *vfa,   /* Dicion�rio anal�tico de flex�es */
  *vfs,   /* Dicion�rio sint�tico de flex�es */
  *vms;   /* Mensagens */

#ifdef MSDOS
/* #include <process.h> */
#define RMODE "rt"
#define WMODE "wt"
#else
#define RMODE "r"
#define WMODE "w"
#endif

/* Internal prototypes: */

int main(int argc, char *argv[]);

FILE *openfile(char *name, char *mode);

bool_t le_linha(byte *lin, int maxlen, long *nlin); 
  /* 
    L� proxima linha de "dfl", terminada por '\n'
    Devolve em "*nlin" o comprimento real da linha 
    (sem contar o '\n' final).  Coloca os bytes lidos
    em *lin (at� o m�ximo de "maxlen" bytes, ignorando o excesso).
    D� erro se a �ltima linha n�o tiver '\n' no fim.
    Resultado � "false" sse "stdin" est� esgotado.
  */
  
void processa_linha(byte *lin, long nlin);
  /* 
    Verifica, analisa, e processa uma linha,
    dados o comprimento total "nlin" e os primeiros MAXLLIN
    caracteres "lin[0..min(nlin,MAXLLIN)-1]".  Ambos excluem
    o NEWLINE final.
  */
  
bool_t testa_delim(byte **r, long *n, char del, byte *lin, long nlin);
  /*
    Verifica se o caracter "**r" � igual a "del"; em caso afirmativo
    incrementa "*r", decrementa "*n", e devolve "true";
    caso contr�rio imprime erro e devolve "false". 
  */
  
void mens_linha(char *msg, long num, byte *lin, long nlin); 
  /*
    Imprime uma mensagem "msg" referente � linha numero "num",
    cujo conteudo (truncado em MAXLLIN bytes) � "*lin" e 
    cujo comprimento � "*nlin".
  */
 
void mens_sintaxe(char *msg, long col, byte *lin, long nlin); 
  /*
    Imprime uma mensagem de erro de sintaxe "msg" referente 
    ao caracter indice "col" da linha numero "num_linhas",
    cujo conteudo (truncado em MAXLLIN bytes) � "*lin" e 
    cujo comprimento � "*nlin".
  */ 
  
void imprime_linha(byte *lin, long nlin);
  /* 
    Imprime a linha "lin" no arquivo "vms", com escapes e elipses,
    para fins de debugagem.
  */
 
void define_tabela(bool_t tb[], char *lets);
  /* 
    Define "tb[c]=true" para todo "c" em "*lets"; demais "false".
  */
  
bool_t cata(
    byte **buf, long *nbuf, 
    byte **res, long *nres,
    bool_t tb[], long maxres,
    byte *lin, long nlin
  );
  /* 
    Retorna em (*res) um apontador para o maior prefixo
    da cadeia (*rbuf)[0..(*nbuf)-1] que � composto de caracteres 
    da tabela "tb"; e devolve em *nres o comprimento desse prefixo.
    Tamb�m atualiza (*rbuf) e (*nbuf) para o restante da cadeia. 
    Se o prefixo tem comprimento maior  que "maxres", imprime
    uma mensagem de erro (com a linha "lin"), incrementa num_erros,
    e devolve "false"; caso contr�rio devolve "true".
  */
  
char *txtcat(char *a, char *b);
  /*
    Devolve uma nova cadeia que cont�m a concatena��o de "a" 
    e "b".
  */
 
#define check(cond,msg) \
  { if(!(cond)) \
      { if (vms) fprintf(vms, "** %s\n", msg); \
        if (vms != stderr) fprintf(stderr, "** %s\n", msg); \
        exit(1); \
      }  \
  }

static bool_t tb_pal[256]; /* Bytes v�lidos em palavras */
static bool_t tb_cat[256]; /* Bytes v�lidos em dados gramaticais */ 
static bool_t tb_ord[256]; /* Bytes v�lidos em d�gitos de ordem */ 

static char digit[36] = "0123456789abcdefghijklmnopqrstuvwxyz";

int main(int argc, char *argv[])
  { byte *lin = (byte*) malloc(MAXLLIN + 1);
    long nlin;
    
    check(argc == 3, "usage: dfltovoc <dflin> <outname>");
    dfl = openfile(argv[1], RMODE);
    vfa = openfile(txtcat(argv[2], ".vfa"), WMODE);
    vfs = openfile(txtcat(argv[2], ".vfs"), WMODE);
    vms = openfile(txtcat(argv[2], ".vms"), WMODE); 
    if (vms == NULL) vms = stderr;
    
    define_tabela(tb_pal, "'abcdefghijklmnopqrstuvwxyz�������������-");
    define_tabela(tb_cat, "ABCDEFGHIJKLMNOPQRSTUVWXYZ?2-.");
    define_tabela(tb_ord, "123456789");
    
    while (le_linha(lin, MAXLLIN, &nlin))
      { if (num_linhas > MAXREG) 
          { fprintf(stderr, "excesso de linhas na entrada\n");
            return(1);
          }
        num_linhas++;
        processa_linha(lin, nlin); 
      }
   
    fprintf(vms, "\n");
    fprintf(vms, "%8ld linhas no dicionario de entrada\n", num_linhas);
    fprintf(vms, "\n");
    fprintf(vms, "%8ld linhas vazias\n", num_vazias);
    fprintf(vms, "%8ld linhas de flex�es inexistentes\n", num_inexistentes);
    fprintf(vms, "%8ld linhas malformadas\n", num_erros);
    fprintf(vms, "\n");
    fprintf(vms, "%8ld palavras v�lidas no vocabul�rio\n", num_palavras);
    fprintf(vms, "%8ld entradas anal�tcas gravadas\n", num_analiticas);
    fprintf(vms, "%8ld entradas sint�ticas gravadas\n", num_sinteticas);
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
  { int c = getc(dfl);
    byte *p = lin;
    long n = 0;
    if (c == EOF)
      { check(ferror(dfl) == 0, "erro de leitura em stdin"); 
        *nlin = 0;
        return(false);
      }
    else
      { while((c != EOF) && (c != '\n'))
          { if (n < maxlen) *p = (byte) c; 
            n++; p++;
            c = getc(dfl);
          }
        if (c == EOF)
          { check(ferror(dfl) == 0, "erro de leitura em stdin"); 
            check(false, "ultima linha de stdin n�o termina com NEWLINE");
          }
        *nlin = n;
        return(true);
      }
  }
  
void mens_linha(char *msg, long num, byte *lin, long nlin)
  { fprintf(vms, "** linha %8ld", num);
    fprintf(vms, ": %s\n", msg);
    if (lin != NULL) imprime_linha(lin, nlin);
    check(num_erros < MAXERR, "excesso de erros");
  }

void processa_linha(byte *lin, long nlin)
  { byte *pal, *ord, *cls, *flx, *can;
    byte *sufpal, *sufcan;
    long npal, nord, ncls, nflx, ncan, nsufpal, nsufcan;
    
    /* Verifica tamanho e caracteres: */
    if (nlin > MAXLLIN) 
      { mens_linha("muito comprida (ignorada)", num_linhas, lin, nlin);
        num_erros++;
        return;
      }
    if (nlin == 0) 
      { mens_linha("vazia (ignorada)", num_linhas, lin, nlin);
        num_vazias++;
        return;
      }

    /* Divide em campos (pal,npal) (cls,ncls),(flx,nflx) (can,ncan) */
    { byte *r = lin; long n = nlin; 
    
      /* Abre par�nteses: */
      if (! testa_delim(&r, &n, '(', lin, nlin)) return;
        
      /* Palavra: */
      if (! cata(&r, &n, &pal, &npal, tb_pal, MAXLPAL, lin, nlin)) return; 
      
      /* Separador: */
      if (! testa_delim(&r, &n, ',', lin, nlin)) return;
      
      /* D�gito de ordem */
      if (! cata(&r, &n, &ord, &nord, tb_ord, MAXLORD, lin, nlin)) return; 
      
      /* Separador: */
      if (! testa_delim(&r, &n, ',', lin, nlin)) return;
      
      /* Classe gramatical */
      if(! cata(&r, &n, &cls, &ncls, tb_cat, MAXLCLS, lin, nlin)) return; 
      
      /* Separador: */
      if (! testa_delim(&r, &n, ',', lin, nlin)) return;

      /* C�digo de flex�o */
      if(! cata(&r, &n, &flx, &nflx, tb_cat, MAXLFLX, lin, nlin)) return; 
      
      /* Separador: */
      if (! testa_delim(&r, &n, ',', lin, nlin)) return;

      /* Forma can�nica: */
      if(! cata(&r, &n, &can, &ncan, tb_pal, MAXLPAL, lin, nlin)) return; 
      
      /* Fecha par�ntese: */
      if (! testa_delim(&r, &n, ')', lin, nlin)) return;

      /* fim de linha: */
      if (n != 0)
        { mens_linha("lixo no fim da linha (ignorada)", num_linhas, lin, nlin);
          num_erros++;
          return;
        }
    }

    /* Campos vazios: */
    if (nord == 0) 
      { mens_linha("ordem vazia (ignorada)", num_linhas, lin, nlin);
        num_erros++;
        return;
      }

    if (ncls == 0) 
      { mens_linha("classe vazia (ignorada)", num_linhas, lin, nlin);
        num_erros++;
        return;
      }

    if (ncan == 0)
      { mens_linha("forma can�nica vazia", num_linhas, lin, nlin);
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
      { mens_linha("sufixo da forma can�nica muito longo", num_linhas, lin, nlin);
        num_erros++;
        return;
      }

    if (npal == 0) 
      { mens_linha("entrada de forma inexistente", num_linhas, lin, nlin);
        num_inexistentes++;
      }
    
    /* Grava resultado: */
    num_palavras++;
    fprintf(vfa, 
      "%.*s%s%.*s%.*s%.*s%s%c%.*s#\n", 
      (int)npal, pal, 
      SEPA1, 
      (int)nord, ord, 
      (int)ncls, cls, 
      (int)nflx, flx, 
      SEPA2, 
      digit[nsufpal], 
      (int)nsufcan, sufcan
    );
    num_analiticas++;
    fprintf(vfs, 
      "%.*s%s%.*s%.*s%s%c%.*s#\n", 
      (int)ncan, can, 
      SEPS1, 
      (int)ncls, cls, 
      (int)nflx, flx, 
      SEPS2, 
      digit[nsufcan], 
      (int)nsufpal, sufpal
    );
    num_sinteticas++;
  } 
  
bool_t testa_delim(byte **buf, long *nbuf, char del, byte *lin, long nlin)
  {
    if (((*nbuf) == 0) || (((char)(**buf)) != del))
      { mens_sintaxe("delimitador n�o encontrado", nlin - (*nbuf), lin, nlin);
        num_erros++;
        return (false);
      } 
    else
      { (*buf)++; (*nbuf)--; return (true); }
  }

bool_t cata(
    byte **buf, long *nbuf, 
    byte **res, long *nres, 
    bool_t tb[], long maxres,
    byte *lin, long nlin
  )
  { (*res) = (*buf);
    (*nres) = 0;
    while (((*nbuf) > 0) && (tb[**buf]))
      { if ((*nres) >= maxres) 
          { mens_sintaxe("campo muito longo (ignorada)", nlin - (*nbuf), lin, nlin);
            num_erros++;
            return (false);
          }
        (*buf)++; (*nbuf)--; (*nres)++; 
      }
    return(true);
  }
    
void mens_sintaxe(char *msg, long col, byte *lin, long nlin)
  { fprintf(vms, "** linha %8ld coluna %3ld", num_linhas, col);
    fprintf(vms, ": %s\n", msg);
    imprime_linha(lin, nlin);
    check(num_erros < MAXERR, "excesso de erros");
  }
  
void imprime_linha(byte *lin, long nlin)
  { long n = nlin; 
    char *ell;
    fprintf(vms, "\"");
    if (n > MAXLLIN) n = MAXLLIN;
    if (n > 40) { n = 40-3; ell = "..."; } else { ell = ""; }
    while (n > 0) 
      { if (((*lin) >= 32) && ((*lin) != 127))
          { fputc((char)(*lin), vms); }
        else
          { fputc((char)(*lin), vms); }
        n--; lin++;
      }
    fprintf(vms, "%s\"\n", ell);
  }
 
void define_tabela(bool_t tb[], char *lets)
  { int k;
    for (k=0; k<256; k++) tb[k] = false;
    while ((*lets)!='\0') {tb[(byte)(*lets)] = true; lets++; }
  }

char *txtcat(char *a, char *b)
  { char *r = malloc(strlen(a) + strlen(b) + 1);
    check(r != NULL, "memory exhausted");
    strcpy(r, a); strcat(r, b); 
    return(r);
  }

  
