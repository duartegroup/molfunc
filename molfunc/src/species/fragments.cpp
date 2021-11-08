#include "iostream"
#include "algorithm"
#include "math.h"
#include "species/fragments.h"
#include "utils.h"


using namespace std;


namespace molfunc{

    Fragment::Fragment() = default;

    Fragment::Fragment(const string& xyz_filename): Molecule(xyz_filename) {
        /*********************************************************
         * Construct a Fragment molecule from a standard .xyz file,
         * which may be of the form:
         *
         *      7
         *      [*]C(C)=O acetyl,come,coch3,ac    <- title
         *      R   1.56360   1.04770 -0.18810
         *      C   0.81620  -0.01900  0.15360
         *      .    .          .          .
         *
         *  where the title line contains the SMILES string and
         *  a set of name aliases of the fragment. R is the atom
         *  that will be deleted in favour of the core molecule.
         *
         *                 O
         *               //
         *          R---C -- Me
         *
         *  Here the index of the dummy (R) atom is 0 and have
         *  aliases acetyl,come,coch3,ac
         *
         * Arguments:
         *      xyz_filename (string):
         ********************************************************/

        // Populate the name aliases of this fragment
        vector<string> smiles_and_aliases = utils::split(xyz_title_line, ' ');

        if (smiles_and_aliases.size() == 2){
            // Assume aliases are second item in the space separated list
            aliases = utils::split(smiles_and_aliases[1], ',');
        }

        if (n_masked_atoms() != 1){
            throw runtime_error("Cannot construct a fragment molecule with "
                                "no or more than one dummy (R) atom");
        }

        this->dummy_idx = masked_atom_idxs()[0];
        this->dummy_nn_idx = graph.first_neighbour(dummy_idx);
    }

    Fragment::Fragment(const Fragment &fragment): Molecule(fragment) {
        // Copy constructor
        this->rot_grid_w = fragment.rot_grid_w;
        this->aliases = fragment.aliases;
        this->dummy_idx = fragment.dummy_idx;
        this->dummy_nn_idx = fragment.dummy_nn_idx;
        this->cached_coordinates = vector<Coordinate>(fragment.coordinates);
    }

    Fragment::Fragment(vector<Atom3D> atoms, vector<string> aliases)
             : Molecule(atoms){
        /************************************************************
         * Construct a fragment from a set of atoms and
         * a list of aliases e.g. Me, CH3...
         ***********************************************************/
        this->aliases = move(aliases);

        if (n_masked_atoms() != 1){
            throw runtime_error("Cannot construct a fragment molecule with "
                                "no or more than one dummy (R) atom");
        }

        this->dummy_idx = masked_atom_idxs()[0];
        this->dummy_nn_idx = graph.first_neighbour(dummy_idx);
    }

    void Fragment::cache_coordinates(){
        /********************************************
         * Cache the coordinates so that they may
         * be reset following a translation/rotation
         *******************************************/

        this->cached_coordinates = vector<Coordinate>(coordinates);
    }

    void Fragment::reset_coordinates(){
        /************************************
         * Reset the coordinates using the
         * cached values
         ***********************************/
        if (cached_coordinates.empty()){
            throw runtime_error("Cannot reset the coordinates, no "
                                "cached coordinates found.");
        }

         for (unsigned long i=0; i<n_atoms(); i++){
             coordinates[i] = cached_coordinates[i];
         }
    }

    void Fragment::rotate(GridPoint &grid_point){
        /************************************************
         * Rotate this fragment using a point in the grid
         ***********************************************/
         rotation_matrix.update(grid_point);
         Species::rotate(rotation_matrix);
    }

    void Fragment::rotate_about_dummy_nn(GridPoint &grid_point){
        /************************************************
         * Rotate this fragment using a point in the grid
         * on a defined atom index (i.e. coordinate in
         * space)
         ***********************************************/
        if (n_unmasked_atoms() == 1) return;

        rotation_matrix.update(grid_point);
        Species::rotate(rotation_matrix, dummy_nn_idx);
    }

    Fragment FragmentLib::fragment(const string& name){
        /*********************************************************
         * Get a fragment from the library given a name, which must
         * match one of the aliases of the fragment
         *
         * Arguments:
         *      name (str):
         *
         * Raises:
         *      (domain_error):
         ********************************************************/
        string l_name = utils::to_lower(name);

        for (auto &fragment : fragments){
            for (auto &alias : fragment.aliases){

                if (alias == l_name){

                    // Need to copy the fragment as it will be modified in-place
                    return Fragment(fragment);
                }
            }// aliases
        }// fragments

        throw domain_error("Failed to find a fragment with an alias "+
                           name+" in the library");
    }

    vector<vector<Fragment>> FragmentLib::fragments_n_repeats(unsigned long n){
        /*********************************************************
         * Generate a
         *
         * Arguments:
         *      n (str): Number of repeats to perform
         *
         * Raises:
         *      (runtime_error):
         ********************************************************/

        auto n_to_generate = (int)pow(FragmentLib::instance().fragments.size(), n);

        if (n_to_generate > 1000){
            throw runtime_error("Tried to generated "+to_string(n_to_generate)+
                                " fragment combinations, unsupported. Must be <1000");
        }

        vector<vector<Fragment>> pools = {};
        pools.reserve(n);

        for (unsigned long i = 0; i < n; i++){
            auto p = fragments;
            pools.push_back(p);
        }

        vector<vector<Fragment>> result = {{}};

        for (auto &pool : pools) {

            if (result.empty()){
                result.push_back(fragments);
                continue;
            }

            vector<vector<Fragment>> tmp = {};
            for (auto y: pool) {
                for (auto x: result) {     // Requires copy?!
                    x.push_back(y);
                    tmp.push_back(x);
                }
            }
            result = tmp;
        }

        return result;
    }
}
/* -----------------------------------------------
 Below data is automatically generated
 -----------------------------------------------*/

namespace molfunc{
     FragmentLib::FragmentLib(){
         fragments = {

         Fragment({
            Atom3D("O", 2.69990, -1.44350, 0.51050),
            Atom3D("C", 2.29350, -0.41620, -0.08530),
            Atom3D("R", 3.10630, 0.41060, -0.73980),
            Atom3D("C", 0.86569, -0.07174, -0.09264),
            Atom3D("C", 0.07802, -0.19636, 1.06213),
            Atom3D("C", -1.28369, 0.12237, 1.01842),
            Atom3D("C", -1.86686, 0.54619, -0.18062),
            Atom3D("C", -1.08641, 0.66327, -1.33401),
            Atom3D("C", 0.27841, 0.36500, -1.28911),
            Atom3D("H", 0.51824, -0.52758, 1.99485),
            Atom3D("H", -1.88878, 0.03685, 1.91144),
            Atom3D("H", -2.92252, 0.78138, -0.21641),
            Atom3D("H", -1.53944, 0.98177, -2.26447),
            Atom3D("H", 0.87162, 0.45071, -2.19174),
         },
         {"benzoyl","coph","coc6h5","bz",})
         ,
         Fragment({
            Atom3D("R", 3.33370, -0.09980, 0.17350),
            Atom3D("O", 2.38300, 0.00780, -0.79880),
            Atom3D("C", 1.04670, 0.01030, -0.40560),
            Atom3D("C", 0.26980, 1.14430, -0.53240),
            Atom3D("C", -1.06040, 1.12380, -0.13380),
            Atom3D("C", -1.66080, -0.00030, 0.39600),
            Atom3D("C", -0.87810, -1.12960, 0.51960),
            Atom3D("C", 0.44380, -1.11730, 0.12600),
            Atom3D("H", 0.71970, 2.05400, -0.95050),
            Atom3D("H", -1.67600, 2.03150, -0.23710),
            Atom3D("H", -2.70370, 0.00670, 0.70130),
            Atom3D("H", -1.31100, -2.03810, 0.93230),
            Atom3D("H", 1.09330, -1.99330, 0.20940),
         },
         {"phenoxy","oc6h5","oph",})
         ,
         Fragment({
            Atom3D("R", 2.51780, -0.66490, -0.26890),
            Atom3D("C", 1.35790, 0.01360, -0.17010),
            Atom3D("O", 0.20060, -0.66890, 0.02950),
            Atom3D("C", -1.03641, -0.06926, 0.13431),
            Atom3D("O", 1.34790, 1.26130, -0.26160),
            Atom3D("H", -1.05159, 0.63001, 0.99765),
            Atom3D("H", -1.80747, -0.85071, 0.29402),
            Atom3D("H", -1.27493, 0.48265, -0.80012),
         },
         {"co2ch3","co2me",})
         ,
         Fragment({
            Atom3D("R", 2.82100, -0.14580, -0.96940),
            Atom3D("C", 2.19060, 0.06540, 0.24740),
            Atom3D("C", 0.70880, 0.03050, 0.08860),
            Atom3D("C", 0.00750, -1.14030, 0.21600),
            Atom3D("C", -1.36770, -1.22000, 0.07690),
            Atom3D("C", -2.09250, -0.06870, -0.20570),
            Atom3D("C", -1.40030, 1.12210, -0.33760),
            Atom3D("C", -0.01460, 1.17210, -0.19220),
            Atom3D("H", 2.49300, 1.03210, 0.69140),
            Atom3D("H", 2.53840, -0.72900, 0.94320),
            Atom3D("H", 0.59380, -2.04080, 0.43920),
            Atom3D("H", -1.90950, -2.16860, 0.18310),
            Atom3D("H", -3.16770, -0.07010, -0.32510),
            Atom3D("H", -1.92950, 2.04950, -0.55920),
            Atom3D("H", 0.52850, 2.11160, -0.29660),
         },
         {"benzyl","ch2ph","ch2c6h5","bn",})
         ,
         Fragment({
            Atom3D("O", -0.15990, 0.40950, -0.00000),
            Atom3D("R", 1.04260, -0.15370, 0.00000),
            Atom3D("H", -0.88270, -0.25570, -0.00000),
         },
         {"hydroxy","oh",})
         ,
         Fragment({
            Atom3D("R", 1.56360, 1.04770, -0.18810),
            Atom3D("C", 0.81620, -0.01900, 0.15360),
            Atom3D("C", -0.63500, 0.01260, -0.11210),
            Atom3D("O", 1.35800, -0.96910, 0.66730),
            Atom3D("H", -0.84390, -0.57920, -1.05060),
            Atom3D("H", -1.20150, -0.51300, 0.71280),
            Atom3D("H", -1.05750, 1.02000, -0.18300),
         },
         {"acetyl","come","coch3","ac",})
         ,
         Fragment({
            Atom3D("O", 1.36520, 1.37560, 0.20820),
            Atom3D("S", 0.86650, -0.04080, 0.12960),
            Atom3D("C", -0.88470, 0.05090, -0.18020),
            Atom3D("F", -1.16740, 0.60680, -1.39810),
            Atom3D("F", -1.41490, 0.81640, 0.86790),
            Atom3D("F", -1.46960, -1.18120, -0.11650),
            Atom3D("R", 1.66670, -0.95140, -0.99100),
            Atom3D("O", 1.03820, -0.67630, 1.48010),
         },
         {"triflate","so2cf3","tf",})
         ,
         Fragment({
            Atom3D("O", 1.23930, 1.41530, -0.30170),
            Atom3D("S", 0.82140, -0.02070, -0.03460),
            Atom3D("C", -0.93950, -0.02770, 0.01760),
            Atom3D("R", 1.50380, -0.52580, 1.41510),
            Atom3D("O", 1.36150, -0.90200, -1.13900),
            Atom3D("H", -1.31800, 0.93100, 0.45140),
            Atom3D("H", -1.30260, -0.03720, -1.03900),
            Atom3D("H", -1.36600, -0.83290, 0.63020),
         },
         {"methanesulfonyl","mesyl","so2me","so2ch3","ms",})
         ,
         Fragment({
            Atom3D("C", -0.06750, 0.00500, -0.01750),
            Atom3D("R", 1.31650, -0.10110, -0.09870),
            Atom3D("H", -0.32500, -0.08740, 1.05350),
            Atom3D("H", -0.37160, 1.00290, -0.37260),
            Atom3D("H", -0.55240, -0.81940, -0.56480),
         },
         {"ch3","methyl","me",})
         ,
         Fragment({
            Atom3D("N", -1.21440, 0.01150, 0.00000),
            Atom3D("C", -0.05840, 0.00260, 0.00000),
            Atom3D("R", 1.27280, -0.01420, 0.00000),
         },
         {"cyano","nitrile","cn",})
         ,
         Fragment({
            Atom3D("R", 3.23910, 0.18580, -1.39480),
            Atom3D("C", 2.08580, -0.30040, -0.86680),
            Atom3D("O", 1.10840, 0.57530, -0.42510),
            Atom3D("C", -0.11300, 0.17990, 0.13360),
            Atom3D("C", -0.88310, 1.44040, 0.47900),
            Atom3D("C", 0.01930, -0.71920, 1.32400),
            Atom3D("C", -0.90940, -0.50630, -0.96920),
            Atom3D("O", 1.83450, -1.52060, -0.74680),
            Atom3D("H", -1.02450, 1.98250, -0.49790),
            Atom3D("H", -0.27850, 2.11170, 1.11110),
            Atom3D("H", -1.88910, 1.23490, 0.86510),
            Atom3D("H", -0.87530, -0.60840, 1.95190),
            Atom3D("H", 0.88630, -0.39230, 1.92970),
            Atom3D("H", 0.08210, -1.75920, 0.98510),
            Atom3D("H", -0.63780, -0.00260, -1.93870),
            Atom3D("H", -2.00190, -0.34200, -0.85090),
            Atom3D("H", -0.64290, -1.55940, -1.08930),
         },
         {"co2tbu","co2(ch3)3","boc",})
         ,
         Fragment({
            Atom3D("Cl", 0.85040, 0.00000, 0.00000),
            Atom3D("R", -0.85040, 0.00000, 0.00000),
         },
         {"chloride","chloro","cl",})
         ,
         Fragment({
            Atom3D("I", 1.01630, 0.00000, 0.00000),
            Atom3D("R", -1.01630, 0.00000, 0.00000),
         },
         {"iodide","iodo","i",})
         ,
         Fragment({
            Atom3D("R", 0.49390, 0.00000, 0.00000),
            Atom3D("H", -0.49390, 0.00000, 0.00000),
         },
         {"hydro","h",})
         ,
         Fragment({
            Atom3D("R", 0.03850, 1.32400, 2.03050),
            Atom3D("O", -0.63620, 0.30740, 1.53310),
            Atom3D("C", -0.17790, -0.01860, 0.25830),
            Atom3D("C", -0.33510, 1.16470, -0.69610),
            Atom3D("C", -0.88030, -1.20140, -0.32040),
            Atom3D("C", 1.31380, -0.28370, 0.34440),
            Atom3D("H", -1.36350, 1.57300, -0.67980),
            Atom3D("H", 0.40020, 1.95310, -0.46500),
            Atom3D("H", -0.05760, 0.76330, -1.70030),
            Atom3D("H", -1.84660, -0.93640, -0.82350),
            Atom3D("H", -1.16380, -1.87260, 0.54220),
            Atom3D("H", -0.22050, -1.77460, -1.00320),
            Atom3D("H", 1.69430, -0.21980, -0.71470),
            Atom3D("H", 1.75850, 0.51930, 0.93310),
            Atom3D("H", 1.47610, -1.29780, 0.76140),
         },
         {"oc3h9","otbu",})
         ,
         Fragment({
            Atom3D("R", 2.68740, -0.39930, 0.02060),
            Atom3D("C", 1.34910, -0.18200, 0.00920),
            Atom3D("C", 0.48930, -1.25340, 0.07990),
            Atom3D("C", -0.88550, -1.10580, 0.07300),
            Atom3D("C", -1.38440, 0.17790, -0.00890),
            Atom3D("C", -0.53650, 1.26490, -0.08050),
            Atom3D("C", 0.84410, 1.10250, -0.07260),
            Atom3D("H", 0.84820, -2.28230, 0.14560),
            Atom3D("H", -1.55810, -1.96120, 0.12940),
            Atom3D("H", -2.44020, 0.37020, -0.01940),
            Atom3D("H", -0.90960, 2.28530, -0.14570),
            Atom3D("H", 1.49610, 1.98320, -0.13060),
         },
         {"phenyl","c6h5","ph",})
         ,
         Fragment({
            Atom3D("R", 1.80220, -0.21200, -0.02340),
            Atom3D("O", 0.82880, 0.66080, 0.05750),
            Atom3D("C", -0.40510, -0.00070, -0.00140),
            Atom3D("H", -0.57360, -0.57170, 0.93340),
            Atom3D("H", -1.17450, 0.81000, -0.10640),
            Atom3D("H", -0.47770, -0.68640, -0.85970),
         },
         {"methoxy","och3","ome",})
         ,
         Fragment({
            Atom3D("C", -1.02560, -1.42460, -0.31970),
            Atom3D("Si", -0.03700, -0.03220, 0.42000),
            Atom3D("C", -0.72300, 1.63210, -0.00110),
            Atom3D("R", -0.11060, -0.21460, 2.11010),
            Atom3D("C", 1.77230, -0.15560, -0.06530),
            Atom3D("H", -0.53760, -1.86620, -1.21750),
            Atom3D("H", -2.06150, -1.14740, -0.55730),
            Atom3D("H", -1.10240, -2.28060, 0.41070),
            Atom3D("H", -1.77790, 1.62220, -0.30480),
            Atom3D("H", -0.16390, 2.14120, -0.82770),
            Atom3D("H", -0.60470, 2.31360, 0.88350),
            Atom3D("H", 1.90150, -1.11110, -0.63290),
            Atom3D("H", 2.03560, 0.70530, -0.72060),
            Atom3D("H", 2.43490, -0.18220, 0.82260),
         },
         {"trimethylsilyl","sime3","si(ch3)3","tms",})
         ,
         Fragment({
            Atom3D("C", 0.76110, -0.01300, -0.00140),
            Atom3D("C", -0.72480, -0.12820, 0.00470),
            Atom3D("R", -1.24770, 1.05610, -0.51870),
            Atom3D("H", 1.06560, 0.46270, 0.96910),
            Atom3D("H", 1.07160, 0.71880, -0.79260),
            Atom3D("H", 1.30330, -0.94430, -0.19140),
            Atom3D("H", -1.09940, -0.15940, 1.06270),
            Atom3D("H", -1.12970, -0.99270, -0.53250),
         },
         {"c2h5","ch2ch3","ethyl","et",})
         ,
         Fragment({
            Atom3D("C", -1.17650, 0.05000, -0.10490),
            Atom3D("N", 0.07490, 0.59830, 0.42310),
            Atom3D("R", 0.08500, 0.25830, 1.74540),
            Atom3D("C", 1.18010, -0.07480, -0.22240),
            Atom3D("H", -1.95220, 0.83050, 0.13680),
            Atom3D("H", -1.16140, -0.02320, -1.20780),
            Atom3D("H", -1.40250, -0.92270, 0.36840),
            Atom3D("H", 0.95970, -1.15440, -0.20780),
            Atom3D("H", 2.13150, 0.11410, 0.31870),
            Atom3D("H", 1.26150, 0.32400, -1.24960),
         },
         {"dimethylamino","nme2",})
         ,
         Fragment({
            Atom3D("N", -0.10560, 0.00950, 0.27720),
            Atom3D("R", 1.19540, -0.06970, -0.07000),
            Atom3D("H", -0.57800, -0.83270, -0.10290),
            Atom3D("H", -0.51180, 0.89280, -0.10430),
         },
         {"amino","nh2",})
         ,
         Fragment({
            Atom3D("C", -0.12260, 1.36940, 0.45770),
            Atom3D("C", 0.02660, 0.08820, -0.29650),
            Atom3D("R", 0.01320, 0.42450, -1.65530),
            Atom3D("C", 1.32110, -0.62360, -0.03970),
            Atom3D("C", -1.15820, -0.84860, -0.09980),
            Atom3D("H", 0.04330, 2.27360, -0.18970),
            Atom3D("H", 0.54010, 1.48390, 1.31720),
            Atom3D("H", -1.16890, 1.49900, 0.84740),
            Atom3D("H", 2.14370, -0.01790, -0.44910),
            Atom3D("H", 1.49340, -0.89480, 1.00490),
            Atom3D("H", 1.27830, -1.56540, -0.65540),
            Atom3D("H", -1.60350, -0.69190, 0.89700),
            Atom3D("H", -1.96600, -0.61680, -0.83700),
            Atom3D("H", -0.84070, -1.87950, -0.30180),
         },
         {"tertbutyl","c(ch3)3","tbu",})
         ,
         Fragment({
            Atom3D("R", 2.98540, 0.64430, -0.41500),
            Atom3D("C", 1.74810, 0.14660, -0.07990),
            Atom3D("O", 0.74450, -0.01050, -1.01610),
            Atom3D("C", -0.43308, -0.69078, -0.78092),
            Atom3D("C", -1.31757, -0.62643, -2.01702),
            Atom3D("O", 1.49080, -0.18390, 1.10140),
            Atom3D("H", -0.20897, -1.75459, -0.54549),
            Atom3D("H", -0.97518, -0.23250, 0.07505),
            Atom3D("H", -2.26038, -1.18264, -1.83247),
            Atom3D("H", -0.79484, -1.08022, -2.88564),
            Atom3D("H", -1.56370, 0.43114, -2.25135),
         },
         {"co2ch2ch3","co2c2h5","co2et",})
         ,
         Fragment({
            Atom3D("Br", 0.92450, 0.00000, 0.00000),
            Atom3D("R", -0.92450, 0.00000, 0.00000),
         },
         {"bromide","bromo","br",})
         ,
         Fragment({
            Atom3D("C", 1.22330, -0.18630, -0.25170),
            Atom3D("C", -0.01330, 0.63390, 0.04160),
            Atom3D("R", -0.35880, 1.26100, -1.15660),
            Atom3D("C", -1.17760, -0.29640, 0.37270),
            Atom3D("H", 1.75780, -0.50940, 0.65280),
            Atom3D("H", 0.85950, -1.08000, -0.83020),
            Atom3D("H", 1.92730, 0.36040, -0.92510),
            Atom3D("H", 0.12240, 1.36750, 0.84050),
            Atom3D("H", -1.09050, -0.72390, 1.38300),
            Atom3D("H", -1.15490, -1.14150, -0.36260),
            Atom3D("H", -2.09530, 0.31490, 0.23570),
         },
         {"isopropyl","ch(ch3)2","ipr",})
         ,
         Fragment({
            Atom3D("R", 0.67290, 0.00000, 0.00000),
            Atom3D("F", -0.67290, 0.00000, 0.00000),
         },
         {"fluoride","fluro","f",})
         ,
         Fragment({
            Atom3D("O", -0.26410, 1.23310, 0.01100),
            Atom3D("N", -0.01890, 0.02170, 0.00030),
            Atom3D("R", -1.02060, -0.85380, -0.00360),
            Atom3D("O", 1.30360, -0.40100, -0.00770),
         },
         {"nitro","no2",})
         ,
         Fragment({
            Atom3D("R", -0.52360, 1.29720, -0.05520),
            Atom3D("C", -0.35490, -0.04930, 0.00590),
            Atom3D("O", 0.94480, -0.48100, -0.08690),
            Atom3D("H", -0.92090, -0.51640, -0.84060),
            Atom3D("H", -0.77490, -0.45420, 0.96120),
            Atom3D("H", 1.62950, 0.20370, 0.01570),
         },
         {"ch2oh",})
         ,
         Fragment({
            Atom3D("C", -1.37990, 2.46080, 0.11870),
            Atom3D("C", -0.67120, 1.17530, 0.05210),
            Atom3D("C", 0.70790, 1.17240, 0.06530),
            Atom3D("C", 1.37520, -0.03140, 0.00300),
            Atom3D("C", 2.86340, 0.01360, 0.02000),
            Atom3D("C", 0.66940, -1.19090, -0.07000),
            Atom3D("C", -0.69000, -1.21850, -0.08480),
            Atom3D("C", -1.42150, -2.49540, -0.16480),
            Atom3D("C", -1.37110, -0.00740, -0.02220),
            Atom3D("R", -2.72680, 0.00570, -0.03460),
            Atom3D("H", -1.87830, 2.72620, -0.82720),
            Atom3D("H", -0.69870, 3.28850, 0.40650),
            Atom3D("H", -2.16700, 2.37300, 0.90120),
            Atom3D("H", 1.24900, 2.11400, 0.12440),
            Atom3D("H", 3.31060, -0.98360, 0.11400),
            Atom3D("H", 3.14290, 0.56970, 0.95070),
            Atom3D("H", 3.24880, 0.61350, -0.82100),
            Atom3D("H", 1.22090, -2.11820, -0.11760),
            Atom3D("H", -1.27170, -3.09720, 0.74620),
            Atom3D("H", -1.01700, -3.09490, -1.02310),
            Atom3D("H", -2.49470, -2.27520, -0.33690),
         },
         {"mesityl","mes",})
         ,
         Fragment({
            Atom3D("R", -1.08540, -0.11590, -0.83880),
            Atom3D("C", -0.00370, 0.01770, -0.00580),
            Atom3D("F", -0.06870, -1.02060, 0.91360),
            Atom3D("F", -0.03180, 1.24810, 0.61440),
            Atom3D("F", 1.18970, -0.12930, -0.68340),
         },
         {"cf3",})
         ,
     };
  }
}
