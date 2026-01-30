collapse_months <- function(months_factor) {
  # Step 1: convert to month numbers
  months_num <- match(as.character(months_factor), month.abb)
  months_num <- sort(unique(months_num))
  
  # EDGE CASE: all months present
  if (length(months_num) == 12) {
    return("Jan – Dec")
  }
  
  # Step 2: detect breaks
  breaks <- c(TRUE, diff(months_num) != 1)  # TRUE marks start of a new sequence
  
  # Step 3: assign group IDs
  group_id <- cumsum(breaks)
  
  # Step 4: split into sequences
  sequences <- split(months_num, group_id)
  
  # Step 5: check wrap-around (last ends at 12 and first starts at 1)
  first_seq <- sequences[[1]]
  last_seq  <- sequences[[length(sequences)]]
  
  if (last_seq[length(last_seq)] == 12 && first_seq[1] == 1) {
    # merge last and first
    sequences[[1]] <- c(last_seq, first_seq)
    sequences <- sequences[-length(sequences)]
  }
  
  # Step 6: convert sequences to strings
  ranges <- sapply(sequences, function(s) {
    if (length(s) == 1) {
      month.abb[s]
    } else {
      paste0(month.abb[s[1]], " – ", month.abb[s[length(s)]])
    }
  })
  
  # Step 7: join with commas
  paste(ranges, collapse = ", ")
}
