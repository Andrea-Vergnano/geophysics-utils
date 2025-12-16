ReadProsys <- function(filename){
# Required dependencies and imports
# install.packages("bitops") # Uncomment to install if needed
library(bitops)         # For bitwise operations
library(utils)          # For file dialog and progress bar
library(tools)          # For file_ext

# Global variables for bit-level reading
bit_buffer <<- 0
bit_buffer_count <<- 0

# Helper function to read unsigned bits from a binary connection
read_ubit <- function(con, n) {
  # Ensure that the bit buffer has at least n bits.
  while(bit_buffer_count < n) {
    # Read one byte (8 bits)
    new_byte_vec <- readBin(con, what = "integer", size = 1, n = 1, signed = FALSE)
    # Check if end of file was reached
    if(length(new_byte_vec) == 0) {
      warning(paste("Reached EOF while trying to read bits for ubit", n))
      # Decide how to handle EOF: return NA, 0, or stop? Let's return NA.
      return(NA)
    }
    new_byte <- new_byte_vec[1]
    # Append new_byte to the right of current bit_buffer
    # Use bitwise shift for potentially better performance and clarity with large numbers
    bit_buffer <<- bitOr(bitShiftL(bit_buffer, 8), new_byte)
    bit_buffer_count <<- bit_buffer_count + 8
  }
  
  # Calculate how many bits to shift right to align the desired bits
  shift_amount <- bit_buffer_count - n
  # Extract the top n bits using right shift
  result <- bitShiftR(bit_buffer, shift_amount)
  
  # Create a mask to remove the top n bits from the buffer
  # mask = (1 << shift_amount) - 1
  mask <- bitFlip(bitShiftL(bitFlip(0), shift_amount)) # Creates a mask of 'shift_amount' ones
  
  # Apply the mask to keep only the lower bits
  bit_buffer <<- bitAnd(bit_buffer, mask)
  bit_buffer_count <<- shift_amount
  
  return(result)
}


# Reset the bit buffer if needed
reset_bit_buffer <- function() {
  bit_buffer <<- 0
  bit_buffer_count <<- 0
}

# %% Written by J. GANCE
# %% V3.0 - 04/11/2019
# %% R code to read IRIS Instruments .bin and .pro files (Translated from MATLAB)
# % Select the .bin or .pro file
# In R we use file.choose for file selection
if(is.na(filename) || nchar(filename) == 0) { # Check if file selection was cancelled or failed
  stop("File selection cancelled or failed.")
} else {
  file <- filename
  fid <- file(file, "rb") # Open in binary read mode
  
  # Determine TailleTotale (total file size)
  file_info <- file.info(filename)
  TailleTotale <- file_info$size
  
  # Move to beginning of file (seek to bof)
  seek(fid, 0, "start")
  
  # Check file extension
  file_ext <- tools::file_ext(filename) # Use tools::file_ext for robustness
  data <- list() # Initialize data list
  
  # --- Read BIN file ---
  if(tolower(file_ext) == "bin") { # Use tolower for case-insensitivity
    message("Reading BIN file...")
    # data.version=fread(fid,1,'ulong'); (ulong is 4 bytes unsigned in MATLAB context here)
    data$version <- readBin(fid, what = "integer", size = 4, n = 1, signed = FALSE)
    # data.TypeOfSyscal=fread(fid,1,'uint8');
    data$TypeOfSyscal <- readBin(fid, what = "integer", size = 1, n = 1, signed = FALSE)
    
    # Read comment block if applicable
    if((data$TypeOfSyscal %in% c(3, 4, 5, 8, 9, 10, 11)) ||
       (data$version >= 2147483650 && (data$TypeOfSyscal %in% c(1, 2, 6, 7)))) {
      # data.comment=fread(fid,1024,'*char');
      # Read raw bytes first to handle potential null terminators correctly
      raw_comment <- readBin(fid, what = "raw", n = 1024)
      null_pos <- which(raw_comment == as.raw(0))
      if (length(null_pos) > 0) {
        raw_comment <- raw_comment[1:(null_pos[1]-1)]
      }
      data$comment <- rawToChar(raw_comment)
    }
    
    # Read ColeCole block if applicable
    if(((data$TypeOfSyscal %in% c(3, 4, 5, 8, 9, 11)) && data$version == 2147483651) ||
       ((data$TypeOfSyscal %in% c(1, 6, 10)) && data$version >= 2147483651)) {
      # data.ColeCole=fread(fid,[64000 3],'float32');
      tmp <- readBin(fid, what = "numeric", size = 4, n = 64000*3, endian = "little")
      if(length(tmp) == 64000*3) { # Check if read was successful
        data$ColeCole <- matrix(tmp, nrow = 64000, ncol = 3, byrow = FALSE)
      } else {
        warning("Could not read expected ColeCole data.")
        data$ColeCole <- matrix(NA_real_, nrow=0, ncol=3) # Assign empty matrix
      }
    }
    
    # Read file paths block if applicable (version >= 0x80000004)
    if(data$version >= 2147483652) {
      # data.CommonFilePath=fread(fid,260,'*char');
      raw_path <- readBin(fid, what = "raw", n = 260)
      null_pos <- which(raw_path == as.raw(0))
      if (length(null_pos) > 0) {
        raw_path <- raw_path[1:(null_pos[1]-1)]
      }
      data$CommonFilePath <- rawToChar(raw_path)
      # data.NbFiles=fread(fid,1,'ushort');
      data$NbFiles <- readBin(fid, what = "integer", size = 2, n = 1, signed = FALSE)
      # Pre-allocate lists for filenames and sizes
      data$SizeFileName <- vector("list", data$NbFiles)
      data$FileNameIabOrVmn <- vector("list", data$NbFiles)
      if (data$NbFiles > 0) { # Check if NbFiles is positive before looping
        for(i_file in 1:data$NbFiles) {
          # data.SizeFileName{i}=fread(fid,1,'ushort');
          size_fn <- readBin(fid, what = "integer", size = 2, n = 1, signed = FALSE)
          data$SizeFileName[[i_file]] <- size_fn
          # data.FileNameIabOrVmn{i}=fread(fid,data.SizeFileName{i},'*char');
          if (size_fn > 0) {
            raw_fn <- readBin(fid, what = "raw", n = size_fn)
            # Assuming filenames might also be null-terminated within their size
            null_pos_fn <- which(raw_fn == as.raw(0))
            if (length(null_pos_fn) > 0) {
              raw_fn <- raw_fn[1:(null_pos_fn[1]-1)]
            }
            data$FileNameIabOrVmn[[i_file]] <- rawToChar(raw_fn)
          } else {
            data$FileNameIabOrVmn[[i_file]] <- "" # Empty filename if size is 0
          }
        }
      }
    }
    
    # Read Measurements block if applicable
    if(data$TypeOfSyscal %in% c(3, 4, 5, 8, 9, 11)) {
      Position <- seek(fid, 0, "current")
      i <- 1 # Measurement counter
      tic <- Sys.time()
      # Initialize progress bar using total file size
      h <- txtProgressBar(min = Position, max = TailleTotale, initial = Position, style = 3)
      
      # Initialize data$Measure list structure
      # Define all potential fields to ensure consistent structure
      data$Measure <- list(
        el_array = list(), MoreTMesure = list(), time = list(), m_dly = list(),
        TypeCpXyz = list(), Q = list(), pos = list(), Ps = list(), Vp = list(),
        In = list(), rho = list(), m = list(), dev = list(), Tm = list(), Mx = list(),
        Channel = list(), NbChannel = list(), Overload = list(), ChannelValide = list(),
        unused = list(), QuadNumber = list(), Name = list(), Latitude = list(),
        Longitude = list(), NbCren = list(), RsChk = list(),
        # Conditional fields initialized as empty lists
        TxVab = list(), TxBat = list(), RxBat = list(), Temperature = list(),
        DateTime = list(), Iabfile = list(), Vmnfile = list()
      )
      
      reset_bit_buffer()  # Reset bit buffer before starting bit-level reading
      
      # While loop to read measurements until EOF
      while(TRUE) { # Loop indefinitely, break condition inside
        current_loop_start_pos <- seek(fid, 0, "current")
        if (current_loop_start_pos >= TailleTotale) break # Check if already at EOF
        
        # Update progress bar
        setTxtProgressBar(h, current_loop_start_pos)
        
        # --- Start reading one measurement record ---
        # Reading electrode array
        el_array_val <- readBin(fid, what = "integer", size = 2, n = 1, signed = TRUE)
        # Check for EOF immediately after first read of the loop
        if(length(el_array_val) == 0) break
        data$Measure$el_array[[i]] <- el_array_val
        
        # Reading ???
        data$Measure$MoreTMesure[[i]] <- readBin(fid, what = "integer", size = 2, n = 1, signed = TRUE)
        # Reading Injection time
        data$Measure$time[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little")
        # Reading M_delay
        data$Measure$m_dly[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little")
        # Reading Kid or not
        typeCpXyz_val <- readBin(fid, what = "integer", size = 2, n = 1, signed = TRUE)
        data$Measure$TypeCpXyz[[i]] <- typeCpXyz_val
        if(!is.na(typeCpXyz_val) && typeCpXyz_val == 0) {
          warning(paste("Measurement", i, ": Data recorded with Syscal KID (TypeCpXyz=0), interpretation might be limited."))
        }
        # Reading ignored parameter
        data$Measure$Q[[i]] <- readBin(fid, what = "integer", size = 2, n = 1, signed = TRUE)
        # Reading electrode positions (12 floats)
        data$Measure$pos[[i]] <- readBin(fid, what = "numeric", size = 4, n = 12, endian = "little")
        # Reading PS
        data$Measure$Ps[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little")
        # Reading Vp
        data$Measure$Vp[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little")
        # Reading In
        data$Measure$In[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little")
        # Reading resistivity
        data$Measure$rho[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little")
        # Reading chargeability
        data$Measure$m[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little")
        # Reading Q (deviation?)
        data$Measure$dev[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little")
        # Reading IP Windows duration (Tm) (20 floats)
        data$Measure$Tm[[i]] <- readBin(fid, what = "numeric", size = 4, n = 20, endian = "little")
        # Reading IP Windows values (Mx) (20 floats)
        data$Measure$Mx[[i]] <- readBin(fid, what = "numeric", size = 4, n = 20, endian = "little")
        
        # Reading bit fields using read_ubit
        data$Measure$Channel[[i]] <- read_ubit(fid, 4)
        data$Measure$NbChannel[[i]] <- read_ubit(fid, 4)
        data$Measure$Overload[[i]] <- read_ubit(fid, 1)
        data$Measure$ChannelValide[[i]] <- read_ubit(fid, 1)
        data$Measure$unused[[i]] <- read_ubit(fid, 6) # Note: Name clashes with PRO file 'unused'
        data$Measure$QuadNumber[[i]] <- read_ubit(fid, 16)
        
        # Reading Name: read 12 raw bytes for null termination handling
        raw_name <- readBin(fid, what = "raw", n = 12)
        null_pos_name <- which(raw_name == as.raw(0))
        if (length(null_pos_name) > 0) {
          raw_name <- raw_name[1:(null_pos_name[1]-1)]
        }
        data$Measure$Name[[i]] <- rawToChar(raw_name)
        
        # Reading Lat/Lon/NbCren/RsChk
        data$Measure$Latitude[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little")
        data$Measure$Longitude[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little")
        data$Measure$NbCren[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little")
        data$Measure$RsChk[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little")
        
        # Conditional reads based on MoreTMesure
        more_t_mesure <- data$Measure$MoreTMesure[[i]] # Get the value for current measurement
        # Initialize conditional fields to NA for this index
        data$Measure$TxVab[[i]] <- NA_real_
        data$Measure$TxBat[[i]] <- NA_real_
        data$Measure$RxBat[[i]] <- NA_real_
        data$Measure$Temperature[[i]] <- NA_real_
        data$Measure$DateTime[[i]] <- NA # Use NA for Date/POSIXct type
        
        if (!is.na(more_t_mesure)) {
          if(more_t_mesure == 2) {
            data$Measure$TxVab[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little")
            data$Measure$TxBat[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little")
            data$Measure$RxBat[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little")
            data$Measure$Temperature[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little")
          } else if(more_t_mesure == 3) {
            data$Measure$TxVab[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little")
            data$Measure$TxBat[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little")
            data$Measure$RxBat[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little")
            data$Measure$Temperature[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little")
            datetime_val <- readBin(fid, what = "numeric", size = 8, n = 1, endian = "little") # Double float
            # Convert Excel/MATLAB date number to R Date
            # Origin for Excel dates is "1899-12-30"
            if (!is.na(datetime_val)) {
              origin_date <- as.Date("1899-12-30")
              # Add the number of days (integer part) to the origin
              # Note: This assumes the value represents days. Time part is ignored.
              data$Measure$DateTime[[i]] <- origin_date + floor(datetime_val)
            }
          }
          # If MoreTMesure is not 2 or 3, fields remain NA
        }
        
        # Conditional reads based on version (for Iabfile/Vmnfile)
        # Initialize to NA first
        data$Measure$Iabfile[[i]] <- NA_integer_
        data$Measure$Vmnfile[[i]] <- NA_integer_
        if(data$version >= 2147483652) { # 0x80000004 en HEXA
          data$Measure$Iabfile[[i]] <- readBin(fid, what = "integer", size = 2, n = 1, signed = TRUE)
          data$Measure$Vmnfile[[i]] <- readBin(fid, what = "integer", size = 2, n = 1, signed = TRUE)
        }
        # --- End reading one measurement record ---
        
        # Check if file position advanced. If not, break to prevent infinite loop.
        current_loop_end_pos <- seek(fid, 0, "current")
        if (current_loop_end_pos == current_loop_start_pos) {
          warning("File position did not advance in BIN reading loop. Breaking.")
          break
        }
        
        i <- i + 1 # Increment measurement counter
      } # End while loop
      
      close(h) # Close progress bar
      toc <- Sys.time()
      message(paste("Finished reading BIN file. Read", i-1, "measurements."))
      print(paste("Elapsed time:", format(toc - tic))) # Use print and format time
      
    } else {
      message("Device type (", data$TypeOfSyscal, ") not managed for BIN measurement reading in this script version.")
    }
    close(fid) # Close file connection
    
    # --- Read PRO file ---
  } else if(tolower(file_ext) == "pro") {
    message("Reading PRO file...")
    # Read header
    data$download_version <- readBin(fid, what = "integer", size = 4, n = 1, signed = TRUE)
    data$size <- readBin(fid, what = "integer", size = 4, n = 1, signed = TRUE)
    page <- -1 # Initialize page variable, used for loop control in original code
    i <- 1 # Measurement counter, start at 1
    
    # Initialize data$Measure list structure for PRO file
    # Define all potential fields read in the loop
    data$Measure <- list(
      page = list(), # Store page number read at start of record
      date = list(), Name = list(), Channel = list(), NbChannel = list(),
      Overload = list(), ChannelValide = list(), unused_bits1 = list(), # Renamed first unused (ubit6)
      QuadNumber = list(), vrunning = list(), vsigned = list(), normalized = list(),
      imperial = list(), unused_bits2 = list(), # Renamed second unused (ubit4)
      timeSet = list(), timeMode = list(), type = list(),
      unused_byte1 = list(), # Renamed third unused (ubit8)
      el_array = list(), time = list(), vdly = list(), mdly = list(),
      tm = list(), # Vector of 20 ushorts
      unused_short1 = list(), # Renamed fourth unused (ushort)
      Latitude = list(), Longitude = list(), inrx = list(), pos = list(), # Vector of 12 floats
      mov = list(), # Vector of 3 floats
      rho = list(), dev = list(), NbCren = list(), RsChk = list(),
      TxVab = list(), TxBat = list(), RxBat = list(), Temperature = list(),
      Ps = list(), Vp = list(), In = list(), m = list(), Mx = list() # Vector of 20 floats
      # rcrc is read but not stored per measurement
    )
    
    reset_bit_buffer()  # Reset bit buffer before starting bit-level reading
    tic <- Sys.time()
    h <- txtProgressBar(min = seek(fid, 0, "current"), max = TailleTotale, style = 3) # Initialize progress bar
    
    repeat { # Loop until break condition (EOF or page == -1)
      current_loop_start_pos <- seek(fid, 0, "current")
      if (current_loop_start_pos >= TailleTotale) break # Check EOF at start
      
      setTxtProgressBar(h, current_loop_start_pos) # Update progress bar
      
      # Read page number
      page_val <- readBin(fid, what = "integer", size = 4, n = 1, signed = TRUE)
      if(length(page_val) == 0) break # Check EOF after reading page
      
      # Check for termination condition from original code
      if(!is.na(page_val) && page_val == -1) {
        message("Encountered page == -1, stopping read as per original logic.")
        break
      }
      data$Measure$page[[i]] <- page_val
      
      # Read date/time stamp (test = fread(fid,1,'uint32');)
      test <- readBin(fid, what = "integer", size = 4, n = 1, signed = FALSE)
      if(length(test) == 0) {
        warning("Reached EOF unexpectedly after reading page number.")
        break # Check for EOF after reading test
      }
      # Interpretation of 'test' (uint32) as date needs clarification. Store as integer for now.
      # Could be seconds since an epoch, or packed date/time.
      data$Measure$date[[i]] <- test
      
      # Read Name (12 chars)
      raw_name_pro <- readBin(fid, what = "raw", n = 12)
      null_pos_name_pro <- which(raw_name_pro == as.raw(0))
      if (length(null_pos_name_pro) > 0) {
        raw_name_pro <- raw_name_pro[1:(null_pos_name_pro[1]-1)]
      }
      data$Measure$Name[[i]] <- rawToChar(raw_name_pro)
      
      # Read bit fields
      data$Measure$Channel[[i]] <- read_ubit(fid, 4)
      data$Measure$NbChannel[[i]] <- read_ubit(fid, 4)
      data$Measure$Overload[[i]] <- read_ubit(fid, 1)
      data$Measure$ChannelValide[[i]] <- read_ubit(fid, 1)
      data$Measure$unused_bits1[[i]] <- read_ubit(fid, 6) # Renamed first unused (ubit6)
      data$Measure$QuadNumber[[i]] <- read_ubit(fid, 16)
      data$Measure$vrunning[[i]] <- read_ubit(fid, 1)
      data$Measure$vsigned[[i]] <- read_ubit(fid, 1)
      data$Measure$normalized[[i]] <- read_ubit(fid, 1)
      data$Measure$imperial[[i]] <- read_ubit(fid, 1)
      data$Measure$unused_bits2[[i]] <- read_ubit(fid, 4) # Renamed second unused (ubit4)
      data$Measure$timeSet[[i]] <- read_ubit(fid, 4)
      data$Measure$timeMode[[i]] <- read_ubit(fid, 4)
      data$Measure$type[[i]] <- read_ubit(fid, 8)
      data$Measure$unused_byte1[[i]] <- read_ubit(fid, 8) # Renamed third unused (ubit8)
      
      # Read remaining fields for the measurement
      data$Measure$el_array[[i]] <- readBin(fid, what = "integer", size = 4, n = 1, signed = TRUE) # int32
      data$Measure$time[[i]] <- readBin(fid, what = "integer", size = 2, n = 1, signed = FALSE) # ushort
      data$Measure$vdly[[i]] <- readBin(fid, what = "integer", size = 2, n = 1, signed = FALSE) # ushort
      data$Measure$mdly[[i]] <- readBin(fid, what = "integer", size = 2, n = 1, signed = FALSE) # ushort
      data$Measure$tm[[i]] <- readBin(fid, what = "integer", size = 2, n = 20, signed = FALSE) # 20 ushorts
      data$Measure$unused_short1[[i]] <- readBin(fid, what = "integer", size = 2, n = 1, signed = FALSE) # Renamed fourth unused (ushort)
      data$Measure$Latitude[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little") # float32
      data$Measure$Longitude[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little") # float32
      data$Measure$inrx[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little") # float32
      data$Measure$pos[[i]] <- readBin(fid, what = "numeric", size = 4, n = 12, endian = "little") # 12 float32
      data$Measure$mov[[i]] <- readBin(fid, what = "numeric", size = 4, n = 3, endian = "little") # 3 float32
      data$Measure$rho[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little") # float32
      data$Measure$dev[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little") # float32
      data$Measure$NbCren[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little") # float32
      data$Measure$RsChk[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little") # float32
      data$Measure$TxVab[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little") # float32
      data$Measure$TxBat[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little") # float32
      data$Measure$RxBat[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little") # float32
      data$Measure$Temperature[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little") # float32
      data$Measure$Ps[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little") # float32
      data$Measure$Vp[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little") # float32
      data$Measure$In[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little") # float32
      data$Measure$m[[i]] <- readBin(fid, what = "numeric", size = 4, n = 1, endian = "little") # float32
      data$Measure$Mx[[i]] <- readBin(fid, what = "numeric", size = 4, n = 20, endian = "little") # 20 float32
      
      # Read record CRC (rcrc = fread(fid,1,'uint32');)
      rcrc <- readBin(fid, what = "integer", size = 4, n = 1, signed = FALSE)
      if(length(rcrc) == 0) {
        warning("Reached EOF while reading rcrc at end of record.")
        break # Check for EOF after reading rcrc
      }
      # rcrc is read but not stored in data$Measure in the original code
      
      # Check if file position advanced, otherwise break to prevent infinite loop
      current_loop_end_pos <- seek(fid, 0, "current")
      if (current_loop_end_pos == current_loop_start_pos) {
        warning("File position did not advance in PRO reading loop. Breaking.")
        break
      }
      
      i <- i + 1 # Increment measurement counter
    } # End repeat loop
    
    close(h) # Close progress bar
    toc <- Sys.time()
    message(paste("Finished reading PRO file. Read", i-1, "measurements."))
    print(paste("Elapsed time:", format(toc - tic)))
    close(fid) # Close file connection
    
  } else {
    # Handle unsupported file types
    close(fid) # Close the file connection opened earlier
    stop(paste("Unsupported file extension:", file_ext, "- Only .bin and .pro are supported."))
  }
  
  # The 'data' list now contains the read information.
  # Measurements are stored in data$Measure, where each field is a list:
  # e.g., data$Measure$el_array[[1]], data$Measure$el_array[[2]], ...
  # For fields reading multiple values per measurement (like pos, Tm, Mx):
  # data$Measure$pos[[1]] is a vector of 12, data$Measure$pos[[2]] is a vector of 12, etc.
  
  # Optional: Convert lists to vectors/matrices for easier analysis in R
  message("Converting measurement lists to vectors/matrices...")
  if (exists("Measure", where = data) && length(data$Measure) > 0 && length(data$Measure[[1]]) > 0) {
     num_measurements <- length(data$Measure[[1]]) # Get number of measurements from first field
     if (num_measurements > 0) {
        for (field_name in names(data$Measure)) {
           # Check if the first element is atomic (not a list/vector itself)
           if (is.atomic(data$Measure[[field_name]][[1]]) && length(data$Measure[[field_name]][[1]]) == 1) {
              data$Measure[[field_name]] <- unlist(data$Measure[[field_name]])
           } else if (is.numeric(data$Measure[[field_name]][[1]])) { # Handle vectors like pos, Tm, Mx
              # Try binding rows; works if all vectors have the same length
              tryCatch({
                 data$Measure[[field_name]] <- do.call(rbind, data$Measure[[field_name]])
              }, error = function(e) {
                 warning(paste("Could not convert field", field_name, "to matrix. Keeping as list."))
              })
           }
           # Add more specific conversions if needed (e.g., for character Name)
        }
     }
  }
  message("Conversion done.")
  
} # End else block for file selection successful

df=as.data.frame(data$Measure)

ERT_raw = data.frame(Xa=df$pos.1,Xb=df$pos.2,Xm=df$pos.3,Xn=df$pos.4,Iab=df$In,Vmn=df$Vp,Rho=df$rho,DevRho=df$dev,Ya=df$pos.5,Yb=df$pos.6,Ym=df$pos.7,Yn=df$pos.8,Za=df$pos.9,Zb=df$pos.10,Zm=df$pos.11,Zn=df$pos.12)

ERT_raw

}
# dati_ert$pseudoSec_x=(dati_ert$Xa+dati_ert$Xb+dati_ert$Xa+dati_ert$Xn)/4
# dati_ert$pseudoSec_z=((dati_ert$Xa+dati_ert$Xm)/2+(dati_ert$Xb+dati_ert$Xn)/2)/2*1/5
# # Clean up global variables if desired
# rm(bit_buffer, bit_buffer_count, envir = .GlobalEnv)
