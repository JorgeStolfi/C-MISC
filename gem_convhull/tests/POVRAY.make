# Makefile para POV-Ray
# Last edited on 2014-07-22 12:40:36 by stolfilocal

# O nome do arquivo principal da descri��o sua cena deve ser 
# OBRIGATORIAMENTE "main.pov", e a imagem final ser� "main.png"

# Os par�metros abaixo definem o tamanho de uma imagem isolada em pixels.
#
FULLWIDTH := 400
FULLHEIGHT := 300

# Os parametros abaixo definem o tamanho da imagem animada:
#
SMALLWIDTH = 40
SMALLHEIGHT = 30

WIDTH := ${FULLWIDTH}
HEIGHT := ${FULLHEIGHT}

# Especifique a densidade de raios por pixel linear. Pode ser NRAYS =
# 1 para rapidez, NRAYS = 2 para qualidade: 
NRAYS := 1

# Liste aqui outros arquivos que fazem parte de seu modelo, como 
# arquivos "#include" chamados pelo seu ".pov",
# imagens usadas em texturas, etc..
#
OTHERINPUTS := 

# Normalmente, voc� n�o deveria precisar mexer nas linhas abaixo:

# Arquivo execut�vel 
POVRAY := povray

# Comando para visualizar as imagens produzidas pelo POV-Ray
IMVIEW := display -title '%d/%f'

# Comando para converter imagens de um formato para outro
CONVERT := convert

# Nome do arquivo ".pov" principal, e prefixo para arqs de sa�da
MAIN := main

# Cuidado: linhas que come�am com 8 colunas em branco devem
# come�ar com 1 TAB, e n�o com SPACEs.

.PHONY: all export

all: ${MAIN}.png
	-${IMVIEW} ${MAIN}.png

${MAIN}.png: ${MAIN}.pov ${OTHERINPUTS}
	-/bin/rm -f ${MAIN}.png
	nice ${POVRAY} \
            +FN +Q9 \
            +W${WIDTH} +H${HEIGHT} \
            +AM1 +A0.0 +R${NRAYS} \
            +D +SP32 +EP4 \
            +L../povinc \
            +L../ttf \
            +l../out \
	    +I${MAIN}.pov \
	    +O${MAIN}.png \
	  3>&2 > ${MAIN}.log

showimage: ${MAIN}.gif
	-${IMVIEW} ${MAIN}.gif
