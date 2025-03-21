# ... existing code ...
combine_csvs_into_excel <- function(csv_dir, excel_file, after_sheet = "Metadata", output_file) {
  
  # Load required package
  if (!requireNamespace("openxlsx2", quietly = TRUE)) {
    stop("Package 'openxlsx2' is required but not installed")
  }
  
  # Get list of CSV files
  csv_files <- list.files(path = csv_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) {
    stop("No CSV files found in directory: ", csv_dir)
  }
  
  # Load Excel workbook
  wb <- openxlsx2::wb_load(excel_file)
  existing_sheets <- openxlsx2::wb_get_sheet_names(wb)
  
  # Handle case sensitivity in sheet name matching
  sheet_index <- which(tolower(existing_sheets) == tolower(after_sheet))
  if (length(sheet_index) == 0) {
    stop("Sheet '", after_sheet, "' not found in workbook (checked case-insensitively)")
  }
  
  # Use the actual case as it appears in the workbook
  after_sheet <- existing_sheets[sheet_index[1]]
  
  # Get index to insert after
  insert_index <- which(existing_sheets == after_sheet)
  
  # Track new sheet names
  new_sheets <- character(0)
  
  # Add each CSV as new sheet
  for (i in seq_along(csv_files)) {
    # Get sheet name from file name with original capitalization preserved
    sheet_name <- tools::file_path_sans_ext(basename(csv_files[i]))
    
    # Check for sheet name length limits (31 characters for Excel)
    if (nchar(sheet_name) > 31) {
      sheet_name <- substr(sheet_name, 1, 31)
      message("Sheet name truncated to 31 characters: ", sheet_name)
    }
    
    # Make sheet name unique if needed
    if (sheet_name %in% c(existing_sheets, new_sheets)) {
      new_name <- paste0(sheet_name, "_", i)
      message("Sheet name '", sheet_name, "' already exists, using '", new_name, "' instead")
      sheet_name <- new_name
    }
    
    # Read and add CSV data
    message("Reading file: ", basename(csv_files[i]))
    data <- utils::read.csv(csv_files[i])
    
    tryCatch({
      # Create the worksheet with explicit name
      message("Adding worksheet: '", sheet_name, "'")
      openxlsx2::wb_add_worksheet(wb, sheet = sheet_name)
      
      # Add data to the worksheet
      message("Adding data to worksheet: '", sheet_name, "'")
      openxlsx2::wb_add_data(wb, sheet = sheet_name, x = data)
      
      # Successfully added, so track the new sheet
      new_sheets <- c(new_sheets, sheet_name)
    }, error = function(e) {
      message("Error adding worksheet '", sheet_name, "': ", e$message)
    })
  }
  
  if (length(new_sheets) == 0) {
    stop("No sheets were successfully added")
  }
  
  # Reorder sheets
  message("Reordering sheets...")
  updated_order <- append(existing_sheets, new_sheets, after = insert_index)
  openxlsx2::wb_set_order(wb, updated_order)
  
  # Add date to output filename
  date_str <- format(Sys.Date(), "%Y%m%d")
  output_path <- paste0(tools::file_path_sans_ext(output_file), 
                        "_", date_str, 
                        ".", tools::file_ext(output_file))
  
  # Save workbook
  message("Saving workbook to: ", output_path)
  openxlsx2::wb_save(wb, output_path, overwrite = TRUE)
  message("CSV files combined and saved to: ", output_path)
}
# ... existing code ...