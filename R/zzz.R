maelstRomboot <- function()
{
  
  OUT <- c(paste0(
    
    r"(                                                                                                ,,ggAAggg,
                                                                                              gAACAAAAAAGAA
                                                                                            /U¨       `*GUCU
                                   _         _      _____                          ,gAUG¨               ¨
                                  | ||      | ||   | ___ \\                       AGUG'   ,gggGGGGGgg
          _ __ ___    __ _   ___  | || ___  | ||_  | |_/ //  ___   _ __ ___      ]GAC   gAAACAAAAUCCCAA
         | '_ ` _ \\ / _` ||/ _ \\| ||/ __||| ___|||  _ <<  / _ \\| '_ ` _ \\     C'  gAAAAU¨'      ¨]C[
         | ||| ||| || ((| || ___//| ||\__ \\| ||_  | | \ \\| (() || ||| ||| ||       CCAGG            *
         |_|||_|||_||\__,_||\___|||_|||___//\____|||_|  \_\\\___//|_|||_|||_||       UCAC    ,aggg,
                                                                                      ]AC   AAAACCAAa
                                                                                            *'-¨*ACAA        
                                                                                                 ]CAC
          G C C                           U G G                           G C A                  CAA
      A U | | | U U       | | |       C C | | | U A       | | |       A G | | | A A            _]U`
| | C | | | | | | | A | | | | | | | G | | | | | | | A | | | | | | | U | | | | | | | C | | |o--´
G C       | | |       C G | | | U A       | | |       C G | | | C A       | | |       A |
                          A C A                           A C A)", 
    "\n",
"\n         Modeller of allele-specific transcriptomics maelstRom, version ", packageVersion("maelstRom"),
"\n",
"\n"))
  
  return(OUT)
}

.onAttach <- function(lib, pkg)
{
  #unlockBinding(".maelstRom", asNamespace("maelstRom")) 
  # startup message
  msg <- maelstRomboot()
  if(!interactive())
    msg[1] <- paste("Package 'maelstRom' version", packageVersion("maelstRom"))
  packageStartupMessage(msg)      
  invisible()
}
