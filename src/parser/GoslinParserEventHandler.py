from .BaseParserEventHandler import BaseParserEventHandler
from domain.LipidAdduct import LipidAdduct
from domain.LipidLevel import LipidLevel
from domain.MolecularFattyAcid import MolecularFattyAcid
from domain.LipidFaBondType import LipidFaBondType
from domain.LipidSpeciesInfo import LipidSpeciesInfo
from domain.LipidSpecies import LipidSpecies
from domain.LipidMolecularSubspecies import LipidMolecularSubspecies
from domain.LipidStructuralSubspecies import LipidStructuralSubspecies
from domain.StructuralFattyAcid import StructuralFattyAcid

class GoslinParserEventHandler(BaseParserEventHandler):
    
    def __init__(self):
        super().__init__()
        self.reset_lipid(None)
        
        self.registered_events["lipid_pre_event"] = self.reset_lipid
        self.registered_events["lipid_post_event"] = self.build_lipid
        
        self.registered_events["hg_cl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_mlcl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_pl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_lpl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_lpl_o_pre_event"] = self.set_head_group_name
        self.registered_events["hg_pl_o_pre_event"] = self.set_head_group_name
        self.registered_events["hg_lsl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_dsl_pre_event"] = self.set_head_group_name
        self.registered_events["ch_pre_event"] = self.set_head_group_name
        self.registered_events["hg_che_pre_event"] = self.set_head_group_name
        self.registered_events["mediator_pre_event"] = self.set_head_group_name
        self.registered_events["hg_mgl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_dgl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_sgl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_tgl_pre_event"] = self.set_head_group_name
        
        self.registered_events["gl_species_pre_event"] = self.set_species_level
        self.registered_events["pl_species_pre_event"] = self.set_species_level
        self.registered_events["chc_pre_event"] = self.set_species_level
        self.registered_events["sl_species_pre_event"] = self.set_species_level
        self.registered_events["mediatorc_pre_event"] = self.set_species_level
        
        self.registered_events["fa2_unsorted_pre_event"] = self.set_molecular_subspecies_level
        self.registered_events["fa3_unsorted_pre_event"] = self.set_molecular_subspecies_level
        self.registered_events["fa4_unsorted_pre_event"] = self.set_molecular_subspecies_level
        
        self.registered_events["lcb_pre_event"] = self.new_lcb
        self.registered_events["fa_pre_event"] = self.new_fa
        self.registered_events["fa_post_event"] = self.append_fa
        
        self.registered_events["ether_pre_event"] = self.add_ether
        self.registered_events["old_hydroxyl_pre_event"] = self.add_old_hydroxyl
        self.registered_events["db_count_pre_event"] = self.add_double_bonds
        self.registered_events["carbon_pre_event"] = self.add_carbon
        self.registered_events["hydroxyl_pre_event"] = self.add_hydroxyl
        
    

    def reset_lipid(self, node):
        self.level = LipidLevel.STRUCTURAL_SUBSPECIES
        self.lipid = None
        self.head_group = ""
        self.lcb = None
        self.fa_list = []
        self.current_fa = None
        

    def set_head_group_name(self, node):
        self.head_group = node.get_text()
        
    
    def set_species_level(self, node):
        self.level = LipidLevel.SPECIES
        
        
    def set_molecular_subspecies_level(self, node):
        self.level = LipidLevel.MOLECULAR_SUBSPECIES
        
        
    def new_fa(self, node):
        if self.level == LipidLevel.SPECIES:
            self.current_fa = LipidSpeciesInfo()
            self.current_fa.level = None
            self.current_fa.num_carbon = None
            self.current_fa.num_hydroxyl = None
            self.current_fa.num_double_bonds = None
            self.current_fa.lipid_FA_bond_type = None
            
        elif self.level == LipidLevel.MOLECULAR_SUBSPECIES:
            self.current_fa = MolecularFattyAcid("FA%i" % (len(self.fa_list) + 1), 2, 0, 0, LipidFaBondType.ESTER, False, -1)
            
        elif self.level == LipidLevel.STRUCTURAL_SUBSPECIES:
            self.current_fa = StructuralFattyAcid("FA%i" % (len(self.fa_list) + 1), 2, 0, 0, LipidFaBondType.ESTER, False, 0)
        
        
    def new_lcb(self, node):
        if self.level == LipidLevel.SPECIES:
            self.lcb = StructuralFattyAcid()
            
        elif self.level == LipidLevel.STRUCTURAL_SUBSPECIES:
            self.lcb = StructuralFattyAcid("LCB", 2, 0, 1, LipidFaBondType.ESTER, True, 1)
            
            
    def append_fa(self, node):
        if self.level == LipidLevel.STRUCTURAL_SUBSPECIES:
            self.current_fa.position = len(self.fa_list) + 1
            
        self.fa_list.append(self.current_fa)
        self.current_fa = None
        
        
    def build_lipid(self, node):
        if self.lcb != None:
            for fa in self.fa_list: fa.position += 1
            self.fa_list = [self.lcb] + self.fa_list
        
        if self.level == LipidLevel.SPECIES:
            self.lipid = LipidSpecies(self.head_group, self.fa_list[0])
            
        elif self.level == LipidLevel.MOLECULAR_SUBSPECIES:
            self.current_fa = LipidMolecularSubspecies(self.head_group, self.fa_list)
            
        elif self.level == LipidLevel.STRUCTURAL_SUBSPECIES:
            self.current_fa = LipidStructuralSubspecies(self.head_group, self.fa_list)
    
        
    def add_ether(self, node):
        ether = node.get_text()
        if ether == "a": self.current_fa.lipid_FA_bond_type = LipidFaBondType.ETHER_PLASMANYL
        elif ether == "p": self.current_fa.lipid_FA_bond_type = LipidFaBondType.ETHER_PLASMANYL
        
        
    def add_old_hydroxyl(self, node):
        old_hydroxyl = node.get_text()
        if old_hydroxyl == "d": self.current_fa.num_hydroxyl = 2
        if old_hydroxyl == "t": self.current_fa.num_hydroxyl = 3
        
        
    def add_double_bonds(self, node):
        self.current_fa.num_double_bonds = node.get_text()
        
        
    def add_carbon(self, node):
        self.current_fa.num_carbon = node.get_text()
        
        
    def add_hydroxyl(self, node):
        self.current_fa.num_hydroxyl = node.get_text()
        
        
        
        
        
        

    """
    public LipidAdduct visitLipid(GoslinParser.LipidContext ctx) {
        GoslinParser.Lipid_eofContext lipid = ctx.lipid_eof();
        Optional<Lipid_pureContext> categoryContext = Optional.ofNullable(lipid.lipid_pure());
        Optional<Adduct_infoContext> adductTermContext = Optional.ofNullable(lipid.adduct_info());

        LipidAdduct la = new LipidAdduct(categoryContext.map((cc) -> {
            return new LipidVisitor().visitLipid_pure(cc);
        }).orElse(LipidSpecies.NONE), adductTermContext.map((t) -> {
            return new AdductVisitor().visitAdduct_info(t);
        }).orElse(Adduct.NONE), "");
        return la;
    }

    private static class LipidVisitor extends GoslinBaseVisitor<LipidSpecies> {

        @Override
        public LipidSpecies visitLipid_pure(GoslinParser.Lipid_pureContext ctx) {
            LipidSpecies lipid = null;
            BitSet bs = new BitSet(5);
            bs.set(LipidCategory.ST.ordinal(), ctx.cholesterol() != null);
            bs.set(LipidCategory.GL.ordinal(), ctx.gl() != null);
            bs.set(LipidCategory.FA.ordinal(), ctx.mediatorc() != null);
            bs.set(LipidCategory.GP.ordinal(), ctx.pl() != null);
            bs.set(LipidCategory.SP.ordinal(), ctx.sl() != null);
            LipidCategory contextCategory = LipidCategory.UNDEFINED;
            switch (bs.cardinality()) {
                case 0:
                    throw new PalinomVisitorException("Parsing context did not contain content for any lipid category. Must contain exactly one of " + Arrays.toString(LipidCategory.values()));
                case 1:
                    contextCategory = LipidCategory.values()[bs.nextSetBit(0)];
                    break;
                default:
                    throw new PalinomVisitorException("Parsing context contained content for more than one lipid category. Must contain exactly one of " + Arrays.toString(LipidCategory.values()));
            }
            switch (contextCategory) {
                case ST:
                    if (ctx.cholesterol().chc() != null) {
                        LipidSpeciesInfo lsi = new LipidSpeciesInfo(LipidLevel.SPECIES, 0, 0, 0, LipidFaBondType.UNDEFINED);
                        lipid = new LipidSpecies(ctx.cholesterol().chc().ch().getText(), LipidCategory.ST, Optional.of(LipidClass.CH), Optional.of(lsi));
                        break;
                    } else if (ctx.cholesterol().che() != null) {
                        lipid = handleChe(ctx.cholesterol().che()).orElse(LipidSpecies.NONE);
                        break;
                    } else {
                        throw new PalinomVisitorException("Unhandled sterol lipid: " + ctx.cholesterol().getText());
                    }
                case GL:
                    lipid = handleGlycerolipid(ctx).orElse(LipidSpecies.NONE);
                    break;
                case FA:
                    lipid = new LipidSpecies(ctx.mediatorc().getText(), LipidCategory.FA, Optional.empty(), Optional.empty());
                    break;
                case GP:
                    lipid = handleGlyceroPhospholipid(ctx).orElse(LipidSpecies.NONE);
                    break;
                case SP:
                    if (ctx.sl().dsl() != null) {
                        lipid = handleDsl(ctx.sl().dsl()).orElse(LipidSpecies.NONE);
                    } else if (ctx.sl().lsl() != null) {
                        lipid = handleLsl(ctx.sl().lsl()).orElse(LipidSpecies.NONE);
                    } else {
                        throw new RuntimeException("Unhandled sphingolipid: " + ctx.sl().getText());
                    }
                    break;
                default:
                    throw new PalinomVisitorException("Unhandled contextCategory: " + contextCategory);
            }
            return lipid;
        }

        private Optional<LipidSpecies> handleGlycerolipid(Lipid_pureContext ctx) throws RuntimeException {
            //glycerophospholipids
            //cardiolipin
            if (ctx.gl().dgl() != null) {
                return handleDgl(ctx.gl().dgl());
            } else if (ctx.gl().mgl() != null) {
                return handleMgl(ctx.gl().mgl());
            } else if (ctx.gl().sgl() != null) {
                return handleSgl(ctx.gl().sgl());
            } else if (ctx.gl().tgl() != null) {
                return handleTgl(ctx.gl().tgl());
            } else {
                throw new PalinomVisitorException("Unhandled context state in GL!");
            }
        }

        private Optional<LipidSpecies> handleDsl(GoslinParser.DslContext dsl) {
            String headGroup = dsl.hg_dslc().getText();
            if (dsl.sl_species() != null) { //species level
                //process species level
                return visitSpeciesLcb(headGroup, dsl.sl_species().lcb());
            } else if (dsl.sl_subspecies() != null) {
                //process subspecies
                if (dsl.sl_subspecies().sorted_fa_separator() != null) {
                    //sorted => StructuralSubspecies
                    return visitStructuralSubspeciesLcb(headGroup, dsl.sl_subspecies().lcb(), Arrays.asList(dsl.sl_subspecies().fa()));
                }
            } else {
                throw new PalinomVisitorException("Unhandled context state in DSL!");
            }
            return Optional.empty();
        }

        private Optional<LipidSpecies> handleLsl(GoslinParser.LslContext lsl) {
            String headGroup = lsl.hg_lslc().getText();
            if (lsl.lcb() != null) { //species / subspecies level
                //process structural sub species level
                return visitStructuralSubspeciesLcb(headGroup, lsl.lcb());
            } else {
                throw new PalinomVisitorException("Unhandled context state in LSL!");
            }
        }

        private static boolean hasElements(List<?> l) {
            return (l == null ? false : !l.isEmpty());
        }

        private Optional<LipidSpecies> handleTgl(GoslinParser.TglContext tgl) {
            String headGroup = tgl.hg_tgl().getText();
            if (tgl.gl_species() != null) { //species level
                //process species level
                return visitSpeciesFas(headGroup, tgl.gl_species().fa());
            } else if (tgl.tgl_subspecies() != null) {
                //process subspecies
                if (tgl.tgl_subspecies().fa3().fa3_sorted() != null) {
                    //sorted => StructuralSubspecies
                    log.info("Building structural subspecies");
                    return visitStructuralSubspeciesFas(headGroup, tgl.tgl_subspecies().fa3().fa3_sorted().fa());
                } else if (tgl.tgl_subspecies().fa3().fa3_unsorted() != null) {
                    //unsorted => MolecularSubspecies
                    log.info("Building molecular subspecies");
                    return visitMolecularSubspeciesFas(headGroup, tgl.tgl_subspecies().fa3().fa3_unsorted().fa());
                }
            } else {
                throw new PalinomVisitorException("Unhandled context state in TGL!");
            }
            return Optional.empty();
        }

        private Optional<LipidSpecies> handleSgl(GoslinParser.SglContext sgl) {
            String headGroup = sgl.hg_sgl().getText();
            if (sgl.gl_species() != null) { //species level
                //process species level
                return visitSpeciesFas(headGroup, sgl.gl_species().fa());
            } else if (sgl.dgl_subspecies() != null) {
                //process subspecies
                if (sgl.dgl_subspecies().fa2().fa2_sorted() != null) {
                    //sorted => StructuralSubspecies
                    return visitStructuralSubspeciesFas(headGroup, sgl.dgl_subspecies().fa2().fa2_sorted().fa());
                } else if (sgl.dgl_subspecies().fa2().fa2_unsorted() != null) {
                    //unsorted => MolecularSubspecies
                    return visitMolecularSubspeciesFas(headGroup, sgl.dgl_subspecies().fa2().fa2_unsorted().fa());
                }
            } else {
                throw new PalinomVisitorException("Unhandled context state in SGL!");
            }
            return Optional.empty();
        }

        private Optional<LipidSpecies> handleChe(GoslinParser.CheContext che) {
            String headGroup = che.hg_chec().getText();
            if (che.fa() != null) {
                return visitStructuralSubspeciesFas(headGroup, Arrays.asList(che.fa()));
            } else {
                throw new PalinomVisitorException("Unhandled context state in ChE!");
            }
        }

        private Optional<LipidSpecies> handleMgl(GoslinParser.MglContext mgl) {
            String headGroup = mgl.hg_mgl().getText();
            if (mgl.fa() != null) {
                return visitStructuralSubspeciesFas(headGroup, Arrays.asList(mgl.fa()));
            } else {
                throw new PalinomVisitorException("Unhandled context state in MGL!");
            }
        }

        private Optional<LipidSpecies> handleDgl(GoslinParser.DglContext dgl) {
            String headGroup = dgl.hg_dgl().getText();
            if (dgl.gl_species() != null) { //species level
                //process species level
                return visitSpeciesFas(headGroup, dgl.gl_species().fa());
            } else if (dgl.dgl_subspecies() != null) {
                //process subspecies
                if (dgl.dgl_subspecies().fa2().fa2_sorted() != null) {
                    //sorted => StructuralSubspecies
                    return visitStructuralSubspeciesFas(headGroup, dgl.dgl_subspecies().fa2().fa2_sorted().fa());
                } else if (dgl.dgl_subspecies().fa2().fa2_unsorted() != null) {
                    //unsorted => MolecularSubspecies
                    return visitMolecularSubspeciesFas(headGroup, dgl.dgl_subspecies().fa2().fa2_unsorted().fa());
                }
            } else {
                throw new PalinomVisitorException("Unhandled context state in DGL!");
            }
            return Optional.empty();
        }

        private Optional<LipidSpecies> handleGlyceroPhospholipid(Lipid_pureContext ctx) throws RuntimeException {
            //glycerophospholipids
            //cardiolipin
            if (ctx.pl().cl() != null) {
                return handleCl(ctx.pl().cl());
            } else if (ctx.pl().dpl() != null) {
                return handleDpl(ctx.pl().dpl());
            } else if (ctx.pl().lpl() != null) {
                return handleLpl(ctx.pl().lpl());
            } else if (ctx.pl().mlcl() != null) {
                return handleMlcl(ctx.pl().mlcl());
            } else if (ctx.pl().pl_o() != null) {
                return handlePlo(ctx.pl().pl_o());
            } else {
                throw new PalinomVisitorException("Unhandled context state in PL!");
            }
        }

        private Optional<LipidSpecies> handlePlo(GoslinParser.Pl_oContext ploc) {
            if (ploc.dpl_o() != null) {
                String headGroup = ploc.dpl_o().hg_pl_oc().getText();
                if (ploc.dpl_o().pl_species() != null) {
                    return visitSpeciesFas(headGroup, ploc.dpl_o().pl_species().fa());
                } else if (ploc.dpl_o().pl_subspecies() != null) {
                    if (ploc.dpl_o().pl_subspecies().fa2().fa2_sorted() != null) {
                        return visitStructuralSubspeciesFas(headGroup, ploc.dpl_o().pl_subspecies().fa2().fa2_sorted().fa());
                    } else if (ploc.dpl_o().pl_subspecies().fa2().fa2_unsorted() != null) {
                        return visitMolecularSubspeciesFas(headGroup, ploc.dpl_o().pl_subspecies().fa2().fa2_unsorted().fa());
                    }
                }
            } else if (ploc.lpl_o() != null) {
                String headGroup = ploc.lpl_o().hg_lpl_oc().getText();
                return visitStructuralSubspeciesFas(headGroup, Arrays.asList(ploc.lpl_o().fa()));
            } else {
                throw new PalinomVisitorException("Unhandled context state in PL O!");
            }
            return Optional.empty();
        }

        private Optional<LipidSpecies> handleCl(GoslinParser.ClContext cl) {
            String headGroup = cl.hg_clc().getText();
            if (cl.pl_species() != null) { //species level
                //process species level
                return visitSpeciesFas(headGroup, cl.pl_species().fa());
            } else if (cl.cl_subspecies() != null) {
                //process subspecies
                if (cl.cl_subspecies().fa4().fa4_sorted() != null) {
                    //sorted => StructuralSubspecies
                    return visitStructuralSubspeciesFas(headGroup, cl.cl_subspecies().fa4().fa4_sorted().fa());
                } else if (cl.cl_subspecies().fa4().fa4_unsorted() != null) {
                    //unsorted => MolecularSubspecies
                    return visitMolecularSubspeciesFas(headGroup, cl.cl_subspecies().fa4().fa4_unsorted().fa());
                }
            } else {
                throw new PalinomVisitorException("Unhandled context state in CL!");
            }
            return Optional.empty();
        }

        private Optional<LipidSpecies> handleMlcl(GoslinParser.MlclContext mlcl) {
            String headGroup = mlcl.hg_mlclc().getText();
            if (mlcl.pl_species() != null) { //species level
                //process species level
                return visitSpeciesFas(headGroup, mlcl.pl_species().fa());
            } else if (mlcl.mlcl_subspecies() != null) {
                //process subspecies
                if (mlcl.mlcl_subspecies().fa3().fa3_sorted() != null) {
                    //sorted => StructuralSubspecies
                    return visitStructuralSubspeciesFas(headGroup, mlcl.mlcl_subspecies().fa3().fa3_sorted().fa());
                } else if (mlcl.mlcl_subspecies().fa3().fa3_unsorted() != null) {
                    //unsorted => MolecularSubspecies
                    return visitMolecularSubspeciesFas(headGroup, mlcl.mlcl_subspecies().fa3().fa3_unsorted().fa());
                }
            } else {
                throw new PalinomVisitorException("Unhandled context state in CL!");
            }
            return Optional.empty();
        }

        private Optional<LipidSpecies> handleDpl(GoslinParser.DplContext dpl) {
            String headGroup = dpl.hg_plc().getText();
            if (dpl.pl_species() != null) { //species level
                //process species level
                return visitSpeciesFas(headGroup, dpl.pl_species().fa());
            } else if (dpl.pl_subspecies() != null) {
                //process subspecies
                if (dpl.pl_subspecies().fa2().fa2_sorted() != null) {
                    //sorted => StructuralSubspecies
                    return visitStructuralSubspeciesFas(headGroup, dpl.pl_subspecies().fa2().fa2_sorted().fa());
                } else if (dpl.pl_subspecies().fa2().fa2_unsorted() != null) {
                    //unsorted => MolecularSubspecies
                    return visitMolecularSubspeciesFas(headGroup, dpl.pl_subspecies().fa2().fa2_unsorted().fa());
                }
            } else {
                throw new PalinomVisitorException("Unhandled context state in PL!");
            }
            return Optional.empty();
        }

        private Optional<LipidSpecies> handleLpl(GoslinParser.LplContext lpl) {
            String headGroup = lpl.hg_lplc().getText();
            //lyso PL has one FA, Species=MolecularSubSpecies=StructuralSubSpecies
            if (lpl.fa() != null) {
                return visitStructuralSubspeciesFas(headGroup, Arrays.asList(lpl.fa()));
            } else {
                throw new PalinomVisitorException("Unhandled context state in PL!");
            }
        }

        private Optional<LipidSpecies> visitSpeciesLcb(String headGroup, GoslinParser.LcbContext lcbContext) {
            return Optional.of(new LipidSpecies(headGroup, getSpeciesInfo(lcbContext)));
        }

        private Optional<LipidSpecies> visitSpeciesFas(String headGroup, GoslinParser.FaContext faContext) {
            return Optional.of(new LipidSpecies(headGroup, getSpeciesInfo(faContext)));
        }

        private Optional<LipidSpecies> visitMolecularSubspeciesFas(String headGroup, List<GoslinParser.FaContext> faContexts) {
            List<MolecularFattyAcid> fas = new LinkedList<>();
            for (int i = 0; i < faContexts.size(); i++) {
                MolecularFattyAcid fa = buildMolecularFa(faContexts.get(i), "FA" + (i + 1));
                fas.add(fa);
            }
            MolecularFattyAcid[] arrs = new MolecularFattyAcid[fas.size()];
            fas.toArray(arrs);
            return Optional.of(new LipidMolecularSubspecies(headGroup, arrs));
        }

        private Optional<LipidSpecies> visitStructuralSubspeciesLcb(String headGroup, GoslinParser.LcbContext lcbContext, List<GoslinParser.FaContext> faContexts) {
            List<StructuralFattyAcid> fas = new LinkedList<>();
            StructuralFattyAcid lcbA = buildStructuralLcb(lcbContext, "LCB", 1);
            fas.add(lcbA);
            for (int i = 0; i < faContexts.size(); i++) {
                StructuralFattyAcid fa = buildStructuralFa(faContexts.get(i), "FA" + (i + 1), i + 2);
                fas.add(fa);
            }
            StructuralFattyAcid[] arrs = new StructuralFattyAcid[fas.size()];
            fas.toArray(arrs);
            return Optional.of(new LipidStructuralSubspecies(headGroup, arrs));
        }

        private Optional<LipidSpecies> visitStructuralSubspeciesFas(String headGroup, List<GoslinParser.FaContext> faContexts) {
            List<StructuralFattyAcid> fas = new LinkedList<>();
            for (int i = 0; i < faContexts.size(); i++) {
                StructuralFattyAcid fa = buildStructuralFa(faContexts.get(i), "FA" + (i + 1), i + 1);
                fas.add(fa);
            }
            StructuralFattyAcid[] arrs = new StructuralFattyAcid[fas.size()];
            fas.toArray(arrs);
            return Optional.of(new LipidStructuralSubspecies(headGroup, arrs));
        }

        private Optional<LipidSpecies> visitStructuralSubspeciesLcb(String headGroup, GoslinParser.LcbContext lcbContext) {
            StructuralFattyAcid fa = buildStructuralLcb(lcbContext, "LCB", 1);
            return Optional.of(new LipidStructuralSubspecies(headGroup, fa));
        }

        private Optional<LipidSpeciesInfo> getSpeciesInfo(GoslinParser.FaContext faContext) {
            if (faContext.fa_pure().heavy() != null) {
                throw new RuntimeException("Heavy label in FA_pure context not implemented yet!");
            }
            LipidFaBondType lfbt = getLipidFaBondType(faContext);
            return Optional.of(new LipidSpeciesInfo(
                    LipidLevel.SPECIES,
                    asInt(faContext.fa_pure().carbon(), 0),
                    asInt(faContext.fa_pure().hydroxyl(), 0),
                    asInt(faContext.fa_pure().db(), 0),
                    lfbt));
        }

        private Optional<LipidSpeciesInfo> getSpeciesInfo(GoslinParser.LcbContext lcbContext) {
            Integer hydroxyl = 0;
            if (lcbContext.old_hydroxyl() != null) {
                switch (lcbContext.old_hydroxyl().getText()) {
                    case "t":
                        hydroxyl = 3;
                        break;
                    case "d":
                        hydroxyl = 2;
                        break;
                    default:
                        throw new PalinomVisitorException("Unsupported old hydroxyl prefix: " + lcbContext.old_hydroxyl().getText());
                }
            } else if (lcbContext.hydroxyl() != null) {
                hydroxyl = asInt(lcbContext.hydroxyl(), 0);
            }
            return Optional.of(new LipidSpeciesInfo(
                    LipidLevel.SPECIES,
                    asInt(lcbContext.carbon(), 0),
                    hydroxyl,
                    asInt(lcbContext.db(), 0), LipidFaBondType.ESTER));
        }
    }

    private static class AdductVisitor extends GoslinBaseVisitor<Adduct> {

        @Override
        public Adduct visitAdduct_info(GoslinParser.Adduct_infoContext ctx) {
            String chargeSign = ctx.charge_sign().getText();
            Integer chargeSignValue = 0;
            switch(chargeSign) {
                case "+":
                    chargeSignValue = 1;
                    break;
                case "-":
                    chargeSignValue = -1;
                    break;
                default:
                    chargeSignValue = 0;
            }
            String adductText = ctx.adduct().getText();
            Adduct adduct = new Adduct("",adductText, Integer.parseInt(ctx.charge().getText()), chargeSignValue);
            return adduct;
        }
    }

    private static LipidFaBondType getLipidFaBondType(GoslinParser.FaContext faContext) throws PalinomVisitorException {
        LipidFaBondType lfbt = LipidFaBondType.ESTER;
        if (faContext.ether() != null) {
            if ("a".equals(faContext.ether().getText())) {
                lfbt = LipidFaBondType.ETHER_PLASMANYL;
            } else if ("p".equals(faContext.ether().getText())) {
                lfbt = LipidFaBondType.ETHER_PLASMENYL;
            } else {
                throw new PalinomVisitorException("Unknown ether context value: " + faContext.ether());
            }
        }
        return lfbt;
    }

    public static MolecularFattyAcid buildMolecularFa(GoslinParser.FaContext ctx, String faName) {
        MolecularFattyAcidBuilder fa = MolecularFattyAcid.molecularFaBuilder();
        LipidFaBondType lfbt = getLipidFaBondType(ctx);
        if (ctx.fa_pure() != null) {
            fa.nCarbon(asInt(ctx.fa_pure().carbon(), 0));
            fa.nHydroxy(asInt(ctx.fa_pure().hydroxyl(), 0));
            if (ctx.fa_pure().db() != null) {
                fa.nDoubleBonds(asInt(ctx.fa_pure().db().db_count(), 0));
                if (ctx.fa_pure().db().db_position() != null) {
                    throw new RuntimeException("Support for double bond positions not implemented yet!");
                }
            }
            fa.lipidFaBondType(lfbt);
            return fa.name(faName).build();

        } else {
            throw new PalinomVisitorException("Uninitialized FaContext!");
        }
    }

    public static <T extends RuleNode> Integer asInt(T context, Integer defaultValue) {
        return maybeMapOr(context, (t) -> {
            return Integer.parseInt(t.getText());
        }, defaultValue);
    }

    public static <T> Optional<T> maybe(T t) {
        return Optional.ofNullable(t);
    }

    public static <T, R> R maybeMapOr(T t, Function<? super T, R> mapper, R r) {
        return maybe(t).map(mapper).orElse(r);
    }

    public static StructuralFattyAcid buildStructuralLcb(GoslinParser.LcbContext ctx, String faName, int position) {
        StructuralFattyAcidBuilder fa = StructuralFattyAcid.structuralFaBuilder();
        fa.nCarbon(asInt(ctx.carbon(), 0));
        fa.nHydroxy(asInt(ctx.hydroxyl(), 0));
        if (ctx.db() != null) {
            fa.nDoubleBonds(asInt(ctx.db().db_count(), 0));
            if (ctx.db().db_position() != null) {
                throw new RuntimeException("Support for double bond positions not implemented yet!");
            }
        }
        fa.lipidFaBondType(LipidFaBondType.ESTER);
        return fa.name(faName).position(position).lcb(true).build();
    }

    public static StructuralFattyAcid buildStructuralFa(GoslinParser.FaContext ctx, String faName, int position) {
        StructuralFattyAcidBuilder fa = StructuralFattyAcid.structuralFaBuilder();
        LipidFaBondType lfbt = getLipidFaBondType(ctx);
        if (ctx.fa_pure() != null) {
            fa.nCarbon(asInt(ctx.fa_pure().carbon(), 0));
            fa.nHydroxy(asInt(ctx.fa_pure().hydroxyl(), 0));
            if (ctx.fa_pure().db() != null) {
                fa.nDoubleBonds(asInt(ctx.fa_pure().db().db_count(), 0));
                if (ctx.fa_pure().db().db_position() != null) {
                    throw new RuntimeException("Support for double bond positions not implemented yet!");
                }
            }
            fa.lipidFaBondType(lfbt);
            return fa.name(faName).position(position).build();

        } else {
            throw new PalinomVisitorException("Uninitialized FaContext!");
        }
    }
    """




    """
    namespace LipidCreator
    {    
        
        [Serializable]
        public class ParserEventHandler : BaseParserEventHandler
        {
            public LipidCreator lipidCreator;
            public Lipid lipid;
            public FattyAcidGroupEnumerator fagEnum;
            public FattyAcidGroup fag;
            public int charge;
            public string adduct;
            public bool sortedSeparator;
            public bool expectsEther;
            public int ethers;
            public bool makeUnsupported = false;
            
            public int heavyIsotope = 0;
            public string heavyElement = "";
            public int heavyCount = 1;
            public ElementDictionary heavyElementCounts = null;
            public ArrayList heavyElementCountList = new ArrayList();
            public string heavyName = "";
            public bool addHeavyPrecursor = false;
        
        
            public ParserEventHandler(LipidCreator _lipidCreator) : base()
            {
                lipidCreator = _lipidCreator;
                resetLipidBuilder(null);
                
                
                
                registeredEvents.Add("lipid_pre_event", resetLipidBuilder);
                registeredEvents.Add("lipid_post_event", lipidPostEvent);
                
                registeredEvents.Add("fa_pre_event", FAPreEvent);
                registeredEvents.Add("fa_post_event", FAPostEvent);
                
                registeredEvents.Add("lcb_pre_event", LCBPreEvent);
                registeredEvents.Add("lcb_post_event", LCBPostEvent);
                registeredEvents.Add("carbon_pre_event", CarbonPreEvent);
                registeredEvents.Add("db_count_pre_event", DB_countPreEvent);
                registeredEvents.Add("hydroxyl_pre_event", HydroxylPreEvent);
                registeredEvents.Add("old_hydroxyl_pre_event", OldHydroxylPreEvent);
                registeredEvents.Add("ether_pre_event", EtherPreEvent);
                
                registeredEvents.Add("gl_pre_event", GLPreEvent);
                registeredEvents.Add("pl_pre_event", PLPreEvent);
                registeredEvents.Add("sl_pre_event", SLPreEvent);
                registeredEvents.Add("cholesterol_pre_event", CholesterolPreEvent);
                registeredEvents.Add("mediator_pre_event", MediatorPreEvent);
                
                registeredEvents.Add("hg_mgl_pre_event", HG_MGLPreEvent);
                registeredEvents.Add("hg_dgl_pre_event", HG_DGLPreEvent);
                registeredEvents.Add("hg_sgl_pre_event", HG_SGLPreEvent);
                registeredEvents.Add("hg_tgl_pre_event", HG_TGLPreEvent);
                
                registeredEvents.Add("hg_cl_pre_event", HG_CLPreEvent);
                registeredEvents.Add("hg_mlcl_pre_event", HG_MLCLPreEvent);
                registeredEvents.Add("hg_pl_pre_event", HG_PLPreEvent);
                registeredEvents.Add("hg_lpl_pre_event", HG_LPLPreEvent);
                registeredEvents.Add("hg_lpl_o_pre_event", HG_LPL_OPreEvent);
                registeredEvents.Add("hg_pl_o_pre_event", HG_PL_OPreEvent);
                
                registeredEvents.Add("hg_lsl_pre_event", HG_LSLPreEvent);
                registeredEvents.Add("hg_dsl_pre_event", HG_DSLPreEvent);
                
                registeredEvents.Add("ch_pre_event", ChPreEvent);
                registeredEvents.Add("hg_che_pre_event", HG_ChEPreEvent);
                
                registeredEvents.Add("dpl_post_event", DPLPostEvent);
                registeredEvents.Add("dpl_o_post_event", DPLPostEvent);
                registeredEvents.Add("sl_post_event", SLPostEvent);
                registeredEvents.Add("mediator_post_event", MediatorPostEvent);
                
                registeredEvents.Add("adduct_pre_event", adductPreEvent);
                registeredEvents.Add("charge_pre_event", chargePreEvent);
                registeredEvents.Add("charge_sign_pre_event", charge_signPreEvent);
                registeredEvents.Add("sorted_fa_separator_pre_event", sortedFASeparatorPreEvent);
                
                registeredEvents.Add("heavy_pre_event", resetHeavy);
                registeredEvents.Add("isotope_pre_event", resetHeavyIsotope);
                registeredEvents.Add("isotope_post_event", addHeavyPrecursorElement);
                registeredEvents.Add("isotope_number_pre_event", addIsotopeNumber);
                registeredEvents.Add("isotope_element_pre_event", addIsotopeElement);
                registeredEvents.Add("isotope_count_pre_event", addIsotopeCount);
                registeredEvents.Add("heavy_hg_post_event", setHGIsotopes);
                registeredEvents.Add("gl_species_pre_event", unsupportedEvent);
                registeredEvents.Add("pl_species_pre_event", unsupportedEvent);
                registeredEvents.Add("sl_species_pre_event", unsupportedEvent);
            }
            
            
            
            public void resetLipidBuilder(Parser.TreeNode node)
            {
                lipid = null;
                fagEnum = null;
                fag = null;
                charge = 0;
                adduct = "";
                sortedSeparator = false;
                expectsEther = false;
                ethers = 0;
                
                heavyIsotope = 0;
                heavyElement = "";
                heavyCount = 1;
                heavyElementCounts = null;
                heavyElementCountList = new ArrayList();
                heavyElementCountList.Add(MS2Fragment.createEmptyElementDict());
                heavyName = "";
                addHeavyPrecursor = false;
            
            }
            
            
            
            
            
            
            
            
            
            public void resetHeavy(Parser.TreeNode node)
            {
                heavyIsotope = 0;
                heavyElement = "";
                heavyCount = 1;
                heavyElementCounts = MS2Fragment.createEmptyElementDict();
                heavyName += node.getText();
                addHeavyPrecursor = true;
            }
            

            
            
            public void resetHeavyIsotope(Parser.TreeNode node)
            {
                heavyIsotope = 0;
                heavyElement = "";
                heavyCount = 1;
            }
            
            
            
            
            public void addIsotopeNumber(Parser.TreeNode node)
            {
                heavyIsotope = Convert.ToInt32(node.getText());
            }
            
            
            
            
            public void addIsotopeElement(Parser.TreeNode node)
            {
                heavyElement = node.getText();
            }
            
            
            
            
            public void setHGIsotopes(Parser.TreeNode node)
            {
                heavyElementCountList[0] = heavyElementCounts;
                heavyName += "HG" + node.getText();
            }
            
            
            
            
            public void addIsotopeCount(Parser.TreeNode node)
            {
                heavyCount = Convert.ToInt32(node.getText());
            }
            
            
            
            
            public void addHeavyPrecursorElement(Parser.TreeNode node)
            {
                if (heavyCount < 1)
                {
                    lipid = null;
                    return;
                }
                
                string key = heavyIsotope.ToString() + heavyElement;
                if (MS2Fragment.ELEMENT_POSITIONS.ContainsKey(key))
                {
                    Molecule m = MS2Fragment.ELEMENT_POSITIONS[key];
                    if (heavyElementCounts.ContainsKey(m))
                    {
                        heavyElementCounts[m] = heavyCount;
                    }
                    else
                    {
                        heavyElementCounts.Add(m, heavyCount);
                    }
                }
                else
                {
                    lipid = null;
                }            
            }
            
            
            
            // handling all events
            public void lipidPostEvent(Parser.TreeNode node)
            {
                // first of all, finish ether PC, PE, LPC, LPE
                // flip fatty acids if necessary
                if (lipid != null && lipid.headGroupNames.Count > 0 && (lipid is Phospholipid))
                {
                    if (!((Phospholipid)lipid).isLyso && !((Phospholipid)lipid).isCL)
                    {
                        bool firstFAHasPlamalogen = false;
                        bool secondFAHasPlamalogen = false;
                        foreach (KeyValuePair<string, bool> kvp in ((Phospholipid)lipid).fag1.faTypes)
                        {
                            firstFAHasPlamalogen |= ((kvp.Key.Equals("FAa") && kvp.Value) || (kvp.Key.Equals("FAp") && kvp.Value));
                        }
                        foreach (KeyValuePair<string, bool> kvp in ((Phospholipid)lipid).fag2.faTypes)
                        {
                            secondFAHasPlamalogen |= ((kvp.Key.Equals("FAa") && kvp.Value) || (kvp.Key.Equals("FAp") && kvp.Value));
                        }
                        
                        // flip fatty acids
                        if (!firstFAHasPlamalogen && secondFAHasPlamalogen)
                        {
                            FattyAcidGroup tmp = ((Phospholipid)lipid).fag1;
                            ((Phospholipid)lipid).fag1 = ((Phospholipid)lipid).fag2;
                            ((Phospholipid)lipid).fag2 = tmp;
                        }
                        
                        else if (firstFAHasPlamalogen && secondFAHasPlamalogen)
                        {
                            lipid = new UnsupportedLipid(lipidCreator);
                        }
                    }
                
                    // check for PE O, PC O, LPE O, LPC O
                    if (lipid != null)
                    {
                        if ((new HashSet<string>(){"PC O", "PE O", "LPC O", "LPE O"}).Contains(lipid.headGroupNames[0]))
                        {
                            string lipidClass = lipid.headGroupNames[0];
                            if (lipidClass.Equals("PC O") || lipidClass.Equals("PE O"))
                            {
                                if (((Phospholipid)lipid).fag1.faTypes["FAp"])
                                {
                                    lipidClass = lipidClass + "-p";
                                }
                                else if (((Phospholipid)lipid).fag1.faTypes["FAa"])
                                {
                                    lipidClass = lipidClass + "-a";
                                }
                            }
                            else if  (lipidClass.Equals("LPC O") || lipidClass.Equals("LPE O"))
                            {
                                if (((Phospholipid)lipid).fag1.faTypes["FAp"])
                                {
                                    lipidClass = lipidClass + "-p";
                                }
                                else if (((Phospholipid)lipid).fag1.faTypes["FAa"])
                                {   
                                    lipidClass = lipidClass + "-a";
                                }
                            }
                            lipid.headGroupNames[0] = lipidClass;
                        }
                        if ((new HashSet<string>(){"PC", "PE", "LPC", "LPE"}).Contains(lipid.headGroupNames[0]) && (((Phospholipid)lipid).fag1.faTypes["FAp"] || ((Phospholipid)lipid).fag1.faTypes["FAa"]))
                        {
                            lipid = null;
                        }
                    }
                }
                
            
            
                if (lipid != null && lipid.headGroupNames.Count > 0 && lipidCreator.headgroups.ContainsKey(lipid.headGroupNames[0]))
                {
                
                    foreach (string adduct in Lipid.ADDUCT_POSITIONS.Keys) lipid.adducts[adduct] = false;
                
                    
                    if (charge != 0)
                    {
                        if (Lipid.ADDUCT_POSITIONS.ContainsKey(adduct) && Lipid.ALL_ADDUCTS[Lipid.ADDUCT_POSITIONS[adduct]].charge == charge && lipidCreator.headgroups[lipid.headGroupNames[0]].adductRestrictions[adduct])
                        {
                            lipid.adducts[adduct] = true;
                        }
                        else
                        {
                            lipid = null;
                        }
                    }
                    else
                    {
                        lipid.adducts[lipidCreator.headgroups[lipid.headGroupNames[0]].defaultAdduct] = true;
                    }
                }
                else 
                {
                    lipid = null;
                }
                
                if (lipid != null && expectsEther && ethers != 1)
                {
                    lipid = null;
                }
                
                if (lipid == null && makeUnsupported)
                {
                    lipid = new UnsupportedLipid(lipidCreator);
                }
                
                
                
                
                if (lipid != null)
                {
                    lipid.onlyHeavyLabeled = 0;
                }
                
                
                
                
                
                // adding heavy labeled isotopes if present
                if (lipid != null && !makeUnsupported && addHeavyPrecursor)
                {
                    ElementDictionary hgDictionary = new ElementDictionary(lipidCreator.headgroups[lipid.headGroupNames[0]].elements);
                    
                    MS2Fragment.updateForHeavyLabeled(hgDictionary, (ElementDictionary)heavyElementCountList[0]);
                    if (!MS2Fragment.validElementDict(hgDictionary)) lipid = null;
                    
                    
                    
                
                    if (lipid != null)
                    {
                        
                        string lipidClass = lipid.headGroupNames[0];
                        heavyElementCountList[0] = hgDictionary;
                        lipid.onlyHeavyLabeled = 1;
                        lipidCreator.addHeavyPrecursor(lipid.headGroupNames[0], heavyName, heavyElementCountList);
                        
                    }
                    
                }
                
            }
            
            
            
            public void unsupportedEvent(Parser.TreeNode node)
            {
                lipid = null;
                makeUnsupported = true;
            }
            
            
            
            
            public void sortedFASeparatorPreEvent(Parser.TreeNode node)
            {
                sortedSeparator = true;
            }
            
            
            
            
            public void GLPreEvent(Parser.TreeNode node)
            {
                lipid = new Glycerolipid(lipidCreator);
                fagEnum = new FattyAcidGroupEnumerator((Glycerolipid)lipid);
            }
            
            
            
            
            public void PLPreEvent(Parser.TreeNode node)
            {
                lipid = new Phospholipid(lipidCreator);
                fagEnum = new FattyAcidGroupEnumerator((Phospholipid)lipid);
            }
            
            
            
            
            public void DPLPostEvent(Parser.TreeNode node)
            {
                if (lipid != null)
                {
                    if (lipid.headGroupNames.Count != 0)
                    {
                        Phospholipid gpl = (Phospholipid)lipid;
                        if (!sortedSeparator && (gpl.fag2.faTypes["FAp"] || gpl.fag2.faTypes["FAa"]))
                        {
                            FattyAcidGroup swapFAG = gpl.fag1;
                            gpl.fag1 = gpl.fag2;
                            gpl.fag2 = swapFAG;
                        }
                        if (gpl.fag2.faTypes["FAp"] || gpl.fag2.faTypes["FAa"])
                        {
                            lipid = null;
                            makeUnsupported = true;
                        }
                    }
                    else
                    {
                        lipid = null;
                    }
                }
            }
            
            
            
            
            public void SLPreEvent(Parser.TreeNode node)
            {
                lipid = new Sphingolipid(lipidCreator);
                ((Sphingolipid)lipid).lcb.hydroxylCounts.Clear();
                ((Sphingolipid)lipid).fag.hydroxylCounts.Clear();
                fagEnum = new FattyAcidGroupEnumerator((Sphingolipid)lipid);
            }
            
            
            
            
            public void CholesterolPreEvent(Parser.TreeNode node)
            {
                lipid = new Cholesterol(lipidCreator);
                fagEnum = new FattyAcidGroupEnumerator((Cholesterol)lipid);
            }
            
            
            
            
            public void MediatorPreEvent(Parser.TreeNode node)
            {
                lipid = new Mediator(lipidCreator);
                string headgroup = node.getText();
                lipid.headGroupNames.Add(headgroup);
            }
            
            
            
            
            public void LCBPreEvent(Parser.TreeNode node)
            {
                fag = ((Sphingolipid)lipid).lcb;
                heavyElementCounts = null;
            }
            
            
            
            
            
            public void LCBPostEvent(Parser.TreeNode node)
            {
                FALCBvalidationCheck();
                if (heavyElementCounts != null)
                {
                    heavyElementCountList.Add(heavyElementCounts);
                }
                else
                {
                    heavyElementCountList.Add(MS2Fragment.createEmptyElementDict());
                }
            }
            
            
            
            
            public void FAPreEvent(Parser.TreeNode node)
            {
                fag = (fagEnum != null && fagEnum.MoveNext()) ? fagEnum.Current : null;
                heavyElementCounts = null;
            }
            
            
            
            
            public void FAPostEvent(Parser.TreeNode node)
            {
                FALCBvalidationCheck();
                if (heavyElementCounts != null)
                {
                    heavyElementCountList.Add(heavyElementCounts);
                }
                else
                {
                    heavyElementCountList.Add(MS2Fragment.createEmptyElementDict());
                }
            }
            
            
            
            
            public void FALCBvalidationCheck()
            {
                // check if created fatty acid is valid
                if (fag != null)
                {
                    if (fag.hydroxylCounts.Count == 0)
                    {
                        fag.hydroxylCounts.Add(0);
                    }
                
                    if (fag.carbonCounts.Count == 0 || fag.doubleBondCounts.Count == 0)
                    {
                        lipid = null;
                    }
                    else if (fag.carbonCounts.Count == 1 && fag.doubleBondCounts.Count == 1)
                    {
                        int carbonLength = (new List<int>(fag.carbonCounts))[0];
                        int doubleBondCount = (new List<int>(fag.doubleBondCounts))[0];
                        
                        int maxDoubleBond = (carbonLength - 1) >> 1;
                        if (doubleBondCount > maxDoubleBond) lipid = null;
                        else if (fag.hydroxylCounts.Count == 1)
                        {
                            int hydroxylCount = (new List<int>(fag.hydroxylCounts))[0];
                            if (carbonLength < hydroxylCount) lipid = null;
                        }
                    }
                    else 
                    {
                        lipid = null;
                    }
                    
                    // check if at least one fatty acid type is enabled
                    int enablesFATypes = 0;
                    foreach(KeyValuePair<string, bool> kvp in fag.faTypes) enablesFATypes += kvp.Value ? 1 : 0;                
                    if (enablesFATypes == 0) lipid = null;
                }
                else 
                {
                    lipid = null;
                }
                
                if (lipid != null && fag != null)
                {
                    // fatty acids with plasmalogens must not contain zero double bonds
                    // it's not a bug, it's a feature ;-)
                    if (fag.faTypes["FAp"] && (new List<int>(fag.doubleBondCounts))[0] == 0)
                    {
                        lipid = null;
                        makeUnsupported = true;
                    }
                }
                
                if (lipid != null && fag != null)
                {
                    foreach(int l in fag.carbonCounts) fag.lengthInfo = Convert.ToString(l);
                    foreach(int db in fag.doubleBondCounts) fag.dbInfo = Convert.ToString(db);
                    foreach(int h in fag.hydroxylCounts) fag.hydroxylInfo = Convert.ToString(h);
                }
            }
            
            public void CarbonPreEvent(Parser.TreeNode node)
            {
                if (fag != null)
                {
                    string carbonCount = node.getText();
                    int carbonCountInt = Convert.ToInt32(carbonCount);
                    if (LipidCreator.MIN_CARBON_LENGTH <= carbonCountInt && carbonCountInt <= LipidCreator.MAX_CARBON_LENGTH) fag.carbonCounts.Add(carbonCountInt);
                    else fag = null;
                }
            }
            
            public void DB_countPreEvent(Parser.TreeNode node)
            {
                if (fag != null)
                {
                    string doubleBondCount = node.getText();
                    int doubleBondCountInt = Convert.ToInt32(doubleBondCount);
                    if (LipidCreator.MIN_DB_LENGTH <= doubleBondCountInt && doubleBondCountInt <= LipidCreator.MAX_DB_LENGTH) fag.doubleBondCounts.Add(doubleBondCountInt);
                    else fag = null;
                }
            }
            
            public void OldHydroxylPreEvent(Parser.TreeNode node)
            {
                if (fag != null)
                {
                    string hydroxylCount = node.getText();
                    if (hydroxylCount == "d") fag.hydroxylCounts.Add(2);
                    else if (hydroxylCount == "t") fag.hydroxylCounts.Add(3);
                    else fag = null;
                }
            }
            
            public void HydroxylPreEvent(Parser.TreeNode node)
            {
                if (fag != null)
                {
                    string hydroxylCount = node.getText();
                    int hydroxylCountInt = Convert.ToInt32(hydroxylCount);
                    if (fag.isLCB && LipidCreator.MIN_LCB_HYDROXY_LENGTH <= hydroxylCountInt && hydroxylCountInt <= LipidCreator.MAX_LCB_HYDROXY_LENGTH) fag.hydroxylCounts.Add(hydroxylCountInt);
                    else if ((lipid is Sphingolipid) && !fag.isLCB && LipidCreator.MIN_SPHINGO_FA_HYDROXY_LENGTH <= hydroxylCountInt && hydroxylCountInt <= LipidCreator.MAX_SPHINGO_FA_HYDROXY_LENGTH) fag.hydroxylCounts.Add(hydroxylCountInt);
                    else if (!(lipid is Sphingolipid) && LipidCreator.MIN_HYDROXY_LENGTH <= hydroxylCountInt && hydroxylCountInt <= LipidCreator.MAX_HYDROXY_LENGTH) fag.hydroxylCounts.Add(hydroxylCountInt);
                    else fag = null;
                }
            }
            
            public void EtherPreEvent(Parser.TreeNode node)
            {
                if (fag != null)
                {
                    List<string> keys = new List<string>(fag.faTypes.Keys);
                    foreach(string faTypeKey in keys) fag.faTypes[faTypeKey] = false;
                
                    string faType = node.getText();
                    fag.faTypes["FA" + faType] = true;
                    ++ethers;
                    /*
                    if ((new HashSet<string>{"LPC O-", "LPE O-", "PC O-", "PE O-"}).Contains(lipid.headGroupNames[0]))
                    {
                        lipid.headGroupNames[0] += faType;
                    }
                    */
                }
            }
            
            public void HG_MGLPreEvent(Parser.TreeNode node)
            {
                if (lipid != null)
                {
                    string headgroup = node.getText();
                    lipid.headGroupNames.Add(headgroup);
                    List<string> keys = new List<string>(((Glycerolipid)lipid).fag2.faTypes.Keys);
                    foreach(string faTypeKey in keys) ((Glycerolipid)lipid).fag2.faTypes[faTypeKey] = false;
                    ((Glycerolipid)lipid).fag2.faTypes["FAx"] = true;
                    keys = new List<string>(((Glycerolipid)lipid).fag3.faTypes.Keys);
                    foreach(string faTypeKey in keys) ((Glycerolipid)lipid).fag3.faTypes[faTypeKey] = false;
                    ((Glycerolipid)lipid).fag3.faTypes["FAx"] = true;
                }
            }
            
            public void HG_DGLPreEvent(Parser.TreeNode node)
            {
                if (lipid != null)
                {
                    string headgroup = node.getText();
                    lipid.headGroupNames.Add(headgroup);
                    List<string> keys = new List<string>(((Glycerolipid)lipid).fag3.faTypes.Keys);
                    foreach(string faTypeKey in keys) ((Glycerolipid)lipid).fag3.faTypes[faTypeKey] = false;
                    ((Glycerolipid)lipid).fag3.faTypes["FAx"] = true;
                }
            }
            
            public void HG_SGLPreEvent(Parser.TreeNode node)
            {
                if (lipid != null)
                {
                    string headgroup = node.getText();
                    lipid.headGroupNames.Add(headgroup);
                    List<string> keys = new List<string>(((Glycerolipid)lipid).fag3.faTypes.Keys);
                    foreach(string faTypeKey in keys) ((Glycerolipid)lipid).fag3.faTypes[faTypeKey] = false;
                    ((Glycerolipid)lipid).fag3.faTypes["FAx"] = true;
                    ((Glycerolipid)lipid).containsSugar = true;
                }
            }
            
            public void HG_TGLPreEvent(Parser.TreeNode node)
            {
                if (lipid != null)
                {
                    string headgroup = node.getText();
                    lipid.headGroupNames.Add(headgroup);
                }
            }
            
            public void HG_CLPreEvent(Parser.TreeNode node)
            {
                if (lipid != null)
                {
                    string headgroup = node.getText();
                    lipid.headGroupNames.Add(headgroup);
                    ((Phospholipid)lipid).isCL = true;
                }
            }
            
            public void HG_MLCLPreEvent(Parser.TreeNode node)
            {
                if (lipid != null)
                {
                    string headgroup = node.getText();
                    lipid.headGroupNames.Add(headgroup);
                    ((Phospholipid)lipid).isCL = true;
                    List<string> keys = new List<string>(((Phospholipid)lipid).fag4.faTypes.Keys);
                    foreach(string faTypeKey in keys) ((Phospholipid)lipid).fag4.faTypes[faTypeKey] = false;
                    ((Phospholipid)lipid).fag4.faTypes["FAx"] = true;
                }
            }
            
            public void HG_PLPreEvent(Parser.TreeNode node)
            {
                if (lipid != null)
                {
                    string headgroup = node.getText();
                    lipid.headGroupNames.Add(headgroup);
                }
            }
            
            public void HG_LPLPreEvent(Parser.TreeNode node)
            {
                if (lipid != null)
                {
                    string headgroup = node.getText();
                    lipid.headGroupNames.Add(headgroup);
                    ((Phospholipid)lipid).isLyso = true;
                }
            }
            
            public void HG_PL_OPreEvent(Parser.TreeNode node)
            {
                if (lipid != null)
                {
                    string headgroup = node.getText() + " O";
                    expectsEther = true;
                    lipid.headGroupNames.Add(headgroup);
                }
            }
            
            public void HG_LPL_OPreEvent(Parser.TreeNode node)
            {
                if (lipid != null)
                {
                    string headgroup = node.getText() + " O";
                    expectsEther = true;
                    lipid.headGroupNames.Add(headgroup);
                    ((Phospholipid)lipid).isLyso = true;
                }
            }
            
            public void HG_LSLPreEvent(Parser.TreeNode node)
            {
                if (lipid != null)
                {
                    string headgroup = node.getText();
                    lipid.headGroupNames.Add(headgroup);
                    ((Sphingolipid)lipid).isLyso = true;
                }
            }
            
            public void HG_DSLPreEvent(Parser.TreeNode node)
            {
                if (lipid != null)
                {
                    string headgroup = node.getText();
                    if (headgroup == "GB3") headgroup = "Hex3Cer";
                    lipid.headGroupNames.Add(headgroup);
                }
            }
            
            public void ChPreEvent(Parser.TreeNode node)
            {
                if (lipid != null)
                {
                    string headgroup = node.getText();
                    lipid.headGroupNames.Add(headgroup);
                    List<string> keys = new List<string>(((Cholesterol)lipid).fag.faTypes.Keys);
                    foreach(string faTypeKey in keys) ((Cholesterol)lipid).fag.faTypes[faTypeKey] = false;
                    ((Cholesterol)lipid).fag.faTypes["FAx"] = true;
                }
            }
            
            
            
            public void HG_ChEPreEvent(Parser.TreeNode node)
            {
                if (lipid != null)
                {
                    string headgroup = node.getText();
                    lipid.headGroupNames.Add(headgroup);
                    ((Cholesterol)lipid).containsEster = true;
                }
            }
            
            
            
            
                
            public void SLPostEvent(Parser.TreeNode node)
            {
                if (lipid != null)
                {
                    if (lipid.headGroupNames.Count == 0)
                    {
                        lipid = null;
                    }
                }
            }
                
            public void MediatorPostEvent(Parser.TreeNode node)
            {
                if (lipid != null)
                {
                    if (lipid.headGroupNames.Count == 0)
                    {
                        lipid = null;
                    }
                }
            }
            
            public void adductPreEvent(Parser.TreeNode node)
            {
                if (lipid != null)
                {
                    adduct = node.getText();
                
                }
            }
            
            public void chargePreEvent(Parser.TreeNode node)
            {
                charge = Convert.ToInt32(node.getText());
            }
            
            public void charge_signPreEvent(Parser.TreeNode node)
            {
                charge *= node.getText() == "-" ? -1 : 1;
            }
        }    
    }
    """