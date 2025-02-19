rproject_setup_script = function(add_folders = T,
                                 custom_folders = NULL,
                                 dl_github_r_template = F,
                                 gitignore_lines_to_add = NULL) {
  ### Setup folder structure
  if (!is.null(custom_folders)) {
    folders = custom_folders
  } else {
    folders = c('data', 'figures', 'scripts', 'reports', 'tables')
  }
  
  if (add_folders) {
    for (folder in folders) {
      if (!file.exists(folder)) {
        dir.create(folder, recursive = TRUE)
      }
    }
  }
  
  ### Setup gitignore
  if (dl_github_r_template) {
    git = read.delim(
      'https://raw.githubusercontent.com/github/gitignore/refs/heads/main/R.gitignore',
      header = F
    )
    
  } else {
    git = read.delim(".gitignore")
  }
  
  if (!is.null(gitignore_lines_to_add)) {
    # Expects a vector
    for (line in gitignore_lines_to_add) {
      git = rbind(git, '')
    }
  }
  
  writeLines(git[, 1], '.gitignore')
}