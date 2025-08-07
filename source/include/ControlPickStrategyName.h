#ifndef CONTROLPICKSTRATEGYNAME_H
#define CONTROLPICKSTRATEGYNAME_H

#include "Enums.h"

// Une stratégie peut être définie par l'utilisateur
// C'est donc soit une valeur énumérée, soit une chaîne de caractère représentant le nom de la classe utilisateur qui peut se trouver dans le JSON
struct ControlPickStrategyName {
private:    
    std::string userStrategyName;
    PredefinedStrategyName name;
public:
    // Pour initialiser des tableaux
    ControlPickStrategyName() = default;
    
    // Méthode pour traduire d'une chaîne de caractère à un ControlPickStrategyName
    // Évite d'exposer les constructeurs qui n'acceptent que des paramètres particuliers
    static ControlPickStrategyName create(const std::string &str);

    PredefinedStrategyName getPredefinedStrategyName() const;
    const std::string &getUserStrategyName() const;    

    bool operator==(PredefinedStrategyName n) const;
private:
    ControlPickStrategyName(PredefinedStrategyName predefinedName);
    ControlPickStrategyName(const std::string &str);
};

std::ostream &operator<<(std::ostream &os, const ControlPickStrategyName &name);

struct ControlPickStrategyNameTranslator {
    using internal_type = std::string;
    using external_type = ControlPickStrategyName;

    boost::optional<external_type> get_value(const internal_type &str);
    boost::optional<internal_type> put_value(external_type value);   
};

namespace boost { namespace property_tree {
        template<>
        struct translator_between<std::string, ControlPickStrategyName> {
            using type = ControlPickStrategyNameTranslator;
        };
    }}

#endif /* CONTROLPICKSTRATEGYNAME_H */
