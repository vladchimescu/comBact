## Output real donor concentrations


data_dir = "~/Documents/embl/screen_K"

don_annot = readxl::read_xlsx( file.path(data_dir, 
                                         "donorconc.xlsx"),
                               sheet = "donors_K", skip = 2)

don_annot = na.omit(don_annot)

drugannot = read.table(file.path(data_dir, "classes_drugs-NT.txt"),
                       sep = "\t",
                       header = T, stringsAsFactors = F)

drugannot = dplyr::rename(drugannot, donor = code)
drugannot = dplyr::select(drugannot, Drug, donor)

donannot = inner_join(don_annot, drugannot)
donannot = dplyr::select(donannot, -Drug) %>%
  dplyr::rename(donor_real_conc = Conc,
                run = Run)
donannot = dplyr::mutate(donannot, donor_real_conc = as.numeric(donor_real_conc))

donconc1 = donannot
donconc2 = dplyr::mutate(donconc1, donor_real_conc = donor_real_conc / 2)
donconc3 = dplyr::mutate(donconc2, donor_real_conc = donor_real_conc / 2)

donconc1$donor_conc = 1
donconc2$donor_conc = 2
donconc3$donor_conc = 3

donconc = bind_rows(donconc1,
                    donconc2,
                    donconc3)

write.table(donconc, file = file.path(data_dir, "donconcK.txt"),
            row.names = F, quote = F, col.names = T,
            sep = "\t")


# for strain D
data_dir = "~/Documents/embl/screen_D"

don_annot = readxl::read_xlsx( file.path(data_dir, 
                                         "donorconc.xlsx"),
                               skip = 2)
don_annot = dplyr::mutate(don_annot, donor_real_conc = as.numeric(donor_real_conc))
don_annot = na.omit(don_annot)
donconc1 = don_annot
donconc2 = dplyr::mutate(donconc1, donor_real_conc = donor_real_conc / 2)
donconc3 = dplyr::mutate(donconc2, donor_real_conc = donor_real_conc / 2)

donconc1$donor_conc = 1
donconc2$donor_conc = 2
donconc3$donor_conc = 3

donconc = bind_rows(donconc1,
                    donconc2,
                    donconc3)

write.table(donconc, file = file.path(data_dir, "donconcD.txt"),
            row.names = F, quote = F, col.names = T,
            sep = "\t")

