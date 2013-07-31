#* Copyright © - 2008-2013 ANR Textométrie - http://textometrie.ens-lyon.fr
#*
#* This file is part of the TXM platform.
#*
#* The TXM platform is free software: you can redistribute it and/or modif y
#* it under the terms of the GNU General Public License as published by
#* the Free Software Foundation, either version 3 of the License, or
#* (at your option) any later version.
#*
#* The TXM platform is distributed in the hope that it will be useful,
#* but WITHOUT ANY WARRANTY; without even the implied warranty of
#* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#* General Public License for more details.
#*
#* You should have received a copy of the GNU General Public License
#* along with the TXM platform.  If not, see <http://www.gnu.org/licenses/>.

## Sylvain Loiseau <sloiseau@u-paris10.fr>
## Lise Vaudor <lise.vaudor@ens-lyon.fr>
## Matthieu Decorde <matthieu.decorde@ens-lyon.fr>

phyper_bis=function(v_a,v_b,v_c,v_d){
  v_s2=rep(NA,length(v_a))
  for (j in 1:length(v_a)){
      a=v_a[j]
      b=v_b[j]
      c=v_c[j]
      d=v_d
      a_tmp=a#+1
      s1=dhyper(a_tmp,b,c,d)
      s_tmp=dhyper(a_tmp+1,b,c,d)
      s2=s1+s_tmp
      a_tmp=a_tmp+1
      while(log(s2)!=log(s1)){
        s1=s2
        a_tmp=a_tmp+1
        s_tmp=dhyper(a_tmp,b,c,d)
        s2=s1+s_tmp
      }
      v_s2[j]=s2
  }
  return(v_s2)
}

# old function format for retro-compatibility
`specificites` <-
  function(lexicaltable, types=NULL, parts=NULL) {
	return(specificities(lexicaltable, types, parts));
}

`specificities` <-
  function(lexicaltable, types=NULL, parts=NULL) {
    spe <- specificities.probabilities(lexicaltable, types, parts);
    spelog <- matrix(0, nrow=nrow(spe), ncol=ncol(spe));
    spelog[spe < 0] <- log10(-spe[spe < 0]);
    spelog[spe > 0] <- abs(log10(spe[spe >0]));
    spelog[spe == 0] <- 0;
    spelog[is.infinite(spe)] <- 0;
    spelog <- round(spelog, digits=4);
    rownames(spelog) <- rownames(spe);
    colnames(spelog) <- colnames(spe);
    #class(spelog) <- "specificities";
    #attr(spelog, "l.t") <- spe;
    return(spelog);
  }

`specificities.probabilities` <-
  function(lexicaltable, types=NULL, parts=NULL) {
    
    # if (!is.numeric(lexicaltable)) stop("The lexical table must contain numeric values.");
    
    rowMargin <- rowSums(lexicaltable); # or "F" (the total frequency of all the types).
    colMargin <- colSums(lexicaltable); # or "T" (the size of the parts).
    F <- sum(lexicaltable);             # The grand total (number of tokens in the corpus).
    
    if (! is.null(types)) {      # Filter on tokens to be considered.
      if(is.character(types)) {  # convert the name of types given with "types" into row index numbers.
        if (is.null(rownames(lexicaltable))) stop("The lexical table has no row names and the \"types\" argument is a character vector.");
        if (! all(types %in% rownames(lexicaltable))) stop(paste(
          "Some requested types are not known in the lexical table: ",
          paste(types[! (types %in% rownames(lexicaltable))], collapse=" "))
        ); 
      } else {
        if (any(types < 1)) stop("The row index must be greater than 0.");
        if (max(types) > nrow(lexicaltable)) stop("Row index must be smaller than the number of rows.");
      }
      lexicaltable <- lexicaltable[types,, drop = FALSE];
      rowMargin <- rowMargin[types];
    }
    
    if (! is.null(parts)) {      # Filter on parts to be considered.
      if(is.character(parts)) {  # convert the name of parts given with "parts" into col index numbers.
        if (is.null(colnames(lexicaltable))) stop("The lexical table has no col names and the \"parts\" argument is a character vector.");
        if (! all(parts %in% colnames(lexicaltable))) stop(paste(
          "Some requested parts are not known in the lexical table: ",
          paste(parts[! (parts %in% colnames(lexicaltable))], collapse=" "))
        ); 
      } else {
        if (max(parts) > ncol(lexicaltable)) stop("Column index must be smaller than the number of cols.");
        if (any(parts < 1)) stop("The col index must be greater than 0.");
      }
      lexicaltable <- lexicaltable[,parts, drop=FALSE];
      colMargin <- colMargin[parts];
    }
    
    if (nrow(lexicaltable) == 0 | ncol(lexicaltable) == 0) {
      stop("The lexical table must contains at least one row and one column.");
    }
    
    specif <- matrix(0.0, nrow=nrow(lexicaltable), ncol=ncol(lexicaltable));
    
    for(i in 1:ncol(lexicaltable)) {    # We proceed the whole lexical table by column (i.e. by part).
      
      whiteDrawn <- lexicaltable[,i];  # The frequencies observed in this part for each type.
      white <- rowMargin;     # The total frequencies in the corpus for each type.
      black <- F-white;       # The total complement frequency in the corpus for each type.
      drawn <- colMargin[i];  # The number of tokens in the part.
      
      
      independance    <- (white * drawn) / F;         # The theoretic frequency of each type.
      specif_negative <- whiteDrawn <  independance;  # index of observed frequencies below the theoretic frequencies.
      specif_positive <- whiteDrawn >= independance;  # index of observed frequencies above the theoretic frequencies.
      
      specif[specif_negative,i] <- -phyper (
        whiteDrawn[specif_negative], white[specif_negative], black[specif_negative], drawn
      );
      
      specif[specif_positive,i] <- phyper_bis (
        whiteDrawn[specif_positive], white[specif_positive], black[specif_positive], drawn
      );
    }

    colnames(specif) <- colnames(lexicaltable);
    rownames(specif) <- rownames(lexicaltable);
    
    return(specif);
  }

`specificities.lexicon` <-
  function(lexicon, sublexicon) {
    spe <- specificities.lexicon.probabilities(lexicon, sublexicon);
    spelog <- matrix(0, nrow=nrow(spe), ncol=ncol(spe));
    spelog[spe < 0] <- log10(-spe[spe < 0]);
    spelog[spe > 0] <- abs(log10(spe[spe >=0]));
    spelog[spe == 0] <- 0;
    spelog[is.infinite(spe)] <- 0;
    spelog <- round(spelog, digits=4);
    rownames(spelog) <- rownames(spe);
    colnames(spelog) <- colnames(spe);
    #class(spelog) <- "specificities";
    #attr(spelog, "l.t") <- spe;
    return(spelog);
  }

`lexiconsToLexicalTable` <- function(lexicon, sublexicon) {
	if (! all(names(sublexicon) %in% names(lexicon))) 
	stop(paste(
          sum(! (names(sublexicon) %in% names(lexicon))),
          "types of the sublexicon not found in the lexicon: ",
          )
      ); 

	sub <- numeric(length(lexicon));
	names(sub) <- names(lexicon);
	sub[names(sublexicon)] <- sublexicon;

	complementary.lexicon <- c(lexicon- sub);
	if (any(complementary.lexicon < 0)) 
		stop("type cannot be more frequent in the sublexicon than in the lexicon");

	lexicaltable <- matrix(c(sub, complementary.lexicon), ncol=2);
	rownames(lexicaltable) <- names(lexicon);
	colnames(lexicaltable) <- c("sublexicon", "complementary");
	return(lexicaltable)
}

`specificities.lexicon.new` <- function(lexicon, sublexicon) {
  lexicaltable <- lexiconsToLexicalTable(lexicon, sublexicon);
  return(specificities(lexicaltable,NULL,NULL));
}

`specificities.lexicon.probabilities` <-
  function(lexicon, sublexicon) {
    if (!is.numeric(lexicon)) stop("The lexicon must contain numeric values.");
    if (!is.numeric(sublexicon)) stop("The sublexicon must contain numeric values.");
    if (is.null(names(lexicon))) stop("The lexicon must contain names.");
    if (is.null(names(sublexicon))) stop("The sub lexicon must contain names.");
    
    if (! all(names(sublexicon) %in% names(lexicon)))
      stop(
        paste(
          "Some requested types of the sublexicon are not known in the lexicon: ",
          paste(names(sublexicon)[! (names(sublexicon) %in% names(lexicon))], collapse=" ")
        )
      ); 
    
    F <- sum(lexicon);
    f <- sum(sublexicon);
    
    # complementary.lexicon <- c(lexicon[names(sublexicon)] - sublexicon, lexicon[!names(lexicon) %in% names(sublexicon)]);
    
    if (F < f) {
      stop("The lexicon cannot be smaller than the sublexicon");
    }
    
    whiteDrawn <- numeric(length(lexicon)); # The frequencies observed in this part for each type.
    names(whiteDrawn) <- names(lexicon);
    whiteDrawn[names(sublexicon)] <- sublexicon;
    white <- lexicon; # The total frequencies in the corpus for each type.
    black <- F-white;          # The total complement frequency in the corpus for each type.
    drawn <- f;     # The number of tokens in the part.
    
    # print(whiteDrawn);
    # print(white);
    # print(black);
    # print(drawn);
    
    independance    <- (white * drawn) / F;         # The theoretic frequency of each type.
    
    specif_negative <- whiteDrawn <  independance;  # index of observed frequencies below the theoretic frequencies.
    specif_positive <- whiteDrawn >= independance;  # index of observed frequencies above the theoretic frequencies.
    
    specif <- double(length(lexicon));
    
    specif[specif_negative] <- phyper (
      whiteDrawn[specif_negative], 
      white[specif_negative], 
      black[specif_negative], 
      drawn
    );
    
    specif[specif_positive] <- phyper_bis (
      whiteDrawn[specif_positive], white[specif_positive], black[specif_positive], drawn
    );
    
    names(specif) <- names(lexicon);
    
    return(specif);
  }

`SpecifTopN` <- function(Specif, N=10, file=NULL) {
	symbol = Specif
	top <- c()
	cols <- colnames(symbol)

	for (i in 1:length(cols)) {
		sorted <- sort(symbol[, cols[i]], decreasing=TRUE, index.return= TRUE)$x
		top <- union(top, names(sorted[1:N]))
		top <- union(top, names(sorted[length(sorted) -N: length(sorted)]))  
	}

	symbol <- symbol[top,]
	if (file != NULL) {
		write.table(symbol, file)
	}
}

#print.specificities(x, line=20, part=1, form=NULL, ...) {
#  if (all(is.null(line, part))) {
#    stop("either a line or a part must be specified");
#  }
#  if (all(!is.null(line, part))) {
#    stop("only a line or a part must be specified");
#  }
#}