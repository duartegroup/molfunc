#ifndef MOLFUNC_ATOMS_H
#define MOLFUNC_ATOMS_H
#include "string"
#include "array"
#include <memory>


using namespace std;

namespace molfunc{

    class Atom {

        private:
            array<double, 3> coord = {0.0, 0.0, 0.0};

        public:
            string symbol = "X";

            bool masked = false;

            void translate(array<double, 3> &vec);

            [[nodiscard]] unsigned int atomic_number() const;
            [[nodiscard]] bool is_dummy() const;

            explicit Atom();
            explicit Atom(const string& symbol,
                          double x,
                          double y,
                          double z);

            double* ptr_x = nullptr;
            double* ptr_y = nullptr;
            double* ptr_z = nullptr;

            double x();
            double y();
            double z();
    };

}

#endif //MOLFUNC_ATOMS_H
