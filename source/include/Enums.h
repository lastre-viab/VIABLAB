#ifndef ENUMS_H
#define ENUMS_H

#include <string>
#include "boost/property_tree/ptree.hpp"

// On suppose value un nombre entier, res une chaîne de caractères
#define ENUM_TO_STRING(ENUM_ID, ...) if (value == ENUM_ID) { res = #ENUM_ID; } else
// On suppose str une chaîne de caractères, res une enum
#define STRING_TO_ENUM(ENUM_ID, ...) if (str == #ENUM_ID) { res = ENUM_ID; } else
// On suppose value un nombre entier, res une enum
#define STRING_VALUE_TO_ENUM(ENUM_ID, ENUM_VALUE) if (str == #ENUM_VALUE) { res = ENUM_ID; } else

#define WRITE_ENUM_VALUE(ENUM_ID, ENUM_VALUE) ENUM_ID = ENUM_VALUE,

#define CONCAT(A, B) A ## B

#define DEFINE_TRANSLATOR(TYPE, VALUES)                                \
		/* Classe servant à traduire une chaîne de caractère */            \
		/* dans le JSON vers une enum */                                   \
		/* On accepte dans le JSON la valeur entière de l'enum (ex:2), */  \
		/* la valeur entière de l'enum sous forme de */                    \
		/* chaîne de caractères (ex:"2") et le nom de l'enum (ex:"VD") */  \
		struct CONCAT(StringEnumTranslator, TYPE) {                        \
		    using internal_type = std::string;                             \
		    using external_type = TYPE;                                    \
		    \
		    boost::optional<TYPE> get_value(const std::string &str) {      \
			boost::optional<TYPE> res;                                 \
			\
			VALUES(STRING_TO_ENUM) {                                   \
			    VALUES(STRING_VALUE_TO_ENUM) {                         \
			    res = boost::none;                                 \
			}                                                      \
			}                                                          \
			return res;                                                \
		    }                                                              \
		    boost::optional<std::string> put_value(TYPE value) {           \
			boost::optional<std::string> res;                          \
			\
			VALUES(ENUM_TO_STRING) {                                   \
			    res = boost::none;                                     \
			}                                                          \
			\
			return res;                                                \
		    }                                                              \
		};                                                                 \
		/* Classe permettant à boost de connaître notre traducteur et */   \
		/* de l'utiliser par défaut quand il doit traduire notre enum */   \
		namespace boost { namespace property_tree {                        \
		template<>                                                 \
		struct translator_between<std::string, TYPE> {             \
		    using type = CONCAT(StringEnumTranslator, TYPE);       \
		};                                                         \
		}}                                                             \


/*
 * Prend en argument un type et les valeurs d'une enumération (au format
 * identique à TYPE_TRAJ_VALUES) et crée l'enumération et les classes
 * associées à la traduction du property tree de boost.
 */
#define DEFINE_ENUM(ENUM_TYPE, VALUES)                                 \
		enum ENUM_TYPE                                                     \
		{                                                                  \
	VALUES(WRITE_ENUM_VALUE)                                       \
		};                                                                 \
		\
		inline std::string toString(ENUM_TYPE value) {                     \
			std::string res;                                               \
			VALUES(ENUM_TO_STRING) {                                       \
				res = "";                                                  \
			}                                                              \
			return res;                                                    \
		}                                                                  \
		\
		DEFINE_TRANSLATOR(ENUM_TYPE, VALUES)
/*!
 * valeurs pour la constante  de  type  de trajectoire.
 *
 * L'ajout d'un argument FUNCTION permet, en utilisant les macros ci-dessus,
 * de minimiser les erreurs et le code écrit lors de l'écriture des
 * conversions de la chaîne caractère ou des entiers vers l'énum et vice-versa.
 *
 * On peut voir qu'en substituant FUNCTION par l'expansion de WRITE_ENUM_VALUES,
 * on obtient une suite de paire VD = 1, VL = 2...ect définissant ainsi l'enum
 *
 * Si on susbstitue par STRING_TO_ENUM, on obtient une suite de
 * if (str == "VD") { res = VD; } else if ... else
 * Si on substitue par ENUM_TO_STRING, on obtient une suite de
 * if (value == VD) { res = "VD";} else if ... else
 * Si on substitue par STRING_VALUE_TO_ENUM, on obtient une suite de
 * if (str == "1" ) { res = VD; } else if ... else
 *
 * Ces 3 dernières macros doivent avoir un cas de gestion d'erreur
 * si la valeur n'est pas dans l'enum.
 *
 * -----------------------------------------------------------------------------
 *
 * Par exemple TYPE_TRAJ_VALUES(STRING_TO_ENUM), avec indentation ajoutée,
 * sera substitué par:
 * if ("VD" == str) {
 *     res = VD;
 * }
 * else if ("VL" == str) {
 *     res = VL;
 * }
 * ...
 * else
 *
 *
 * Le "else" devant être complété, l'utilisation correcte de la macro est donc:
 * TYPE_TRAJ_ENUM(STRING_TO_ENUM) {
 *    // Gestion de l'erreur
 * }
 */
#define TYPE_TRAJ_VALUES(FUNCTION)                                        \
		FUNCTION(VD, 1)    /* viable par defaut */                            \
		FUNCTION(VL, 2)    /* viable lourd */                                 \
		FUNCTION(OP, 3)    /* optimale */                                     \
		FUNCTION(VMM, 4)   /* viable micro-macro, minimisant la valeur */     \
		FUNCTION(VDI, 5)   /* viable diversifiant le controle */              \
		FUNCTION(VG, 6)    /* viable garanti */                               \
		FUNCTION(STOCHASTIC, 7)   /* choix des contrôles stochastique */      \
		FUNCTION(WEIGHTED_CONTROLS, 8)  /* choix des contrôles avec pondération */ \
		FUNCTION(CAUTIOUS, 9) /* trajectoire avec bulle évitant les bords*/   \
		FUNCTION(WEIGHTED_CONTROLS_CAUTIOUS, 10) /* CAUTIOUS + WEIGHTED_CONTROLS */\
		FUNCTION(CAUTIOUS_HEAVY, 11) /* CAUTIOUS mais avec un contrôle lourd */ \
		FUNCTION(WEIGHTED_CONTROLS_CAUTIOUS_HEAVY, 12) /* CAUTIOUS_HEAVY + WEIGHTED_CONTROLS */ \
		FUNCTION(STRATEGY_LIST, 13) /* L'algorithme est une liste de stratégies appliquées à chaque prériode choix de contrôle */ \

DEFINE_ENUM(TypeTraj, TYPE_TRAJ_VALUES)

inline bool isViableDefault(TypeTraj type)
    {
    return type == VD || type == STOCHASTIC || type == WEIGHTED_CONTROLS;
    }

inline bool isCautious(TypeTraj type)
    {
    return type == CAUTIOUS || type == WEIGHTED_CONTROLS_CAUTIOUS || type == CAUTIOUS_HEAVY || type == WEIGHTED_CONTROLS_CAUTIOUS_HEAVY;
    }

#define TIME_DISCRETIZATION_SCHEME_VALUES(FUNCTION) \
		FUNCTION(NO_DISCRETIZATION_SCHEME, 0)           \
		FUNCTION(EL, 1)                                 \
		FUNCTION(RK2, 2)                                \
		FUNCTION(RK4, 3)

DEFINE_ENUM(TimeDiscretizationScheme, TIME_DISCRETIZATION_SCHEME_VALUES)

/*!
 * valeurs pour la constante  de  type  de probleme
 */
#define SET_TYPE_VALUES(FUNCTION) \
		FUNCTION(VIAB, 1)             \
		FUNCTION(CAPT, 2)             \
		FUNCTION(VIABG, 3)

DEFINE_ENUM(SetType, SET_TYPE_VALUES)

/*
 * Valeurs pour la constante qui définit la méthode de
 * représentation de l'ensemble
 */
#define GRID_METHOD_VALUES(FUNCTION)                     \
		FUNCTION(BS, 1) /* BitSet */             \
		FUNCTION(MM, 2) /* MicroMacro */         \
		FUNCTION(HBS, 3) /* Hybrid BitSet */     \
		FUNCTION(HMM, 4) /* Hybrid MicroMacro */

DEFINE_ENUM(GridMethod, GRID_METHOD_VALUES)

/*!
 * valeurs pour la constante  qui definit le type de dynamique
 */
#define DYN_TYPE_VALUES(FUNCTION)                                  \
		FUNCTION(CC, 1) /*! continuous time and space dynamics */      \
		FUNCTION(DC, 2) /*! discrete time continuous space dynamics*/  \
		FUNCTION(DD, 3) /*! full discrete ( time and space) dynamics*/ \
		FUNCTION(DH, 4) /*! discrete time hybrid space  dynamics*/     \
		FUNCTION(CH, 5) /*! continous time hybrid space dynamics*/

DEFINE_ENUM(DynType, DYN_TYPE_VALUES)

/*
 *  valeurs pour les constantes définissant la méthode de calcul
 *  de la constante de Lipschitz et de la constante M de borne de la dynamique
 */
#define COMPUTE_METHOD_VALUES(FUNCTION) \
		FUNCTION(ANALYTICAL, 0)             \
		FUNCTION(ANALYTICAL_CALC, 1)        \
		FUNCTION(NUMERICAL_CALC, 2)

DEFINE_ENUM(ComputeMethod, COMPUTE_METHOD_VALUES)

#define FD_DYN_TYPE_VALUES(FUNCTION) \
		FUNCTION(FUNC, 1)                \
		FUNCTION(RETRO, 2)

DEFINE_ENUM(FdDynType, FD_DYN_TYPE_VALUES)

#define TARGET_OR_DEPARTURE_VALUES(FUNCTION) \
		FUNCTION(TARGET, 0)                      \
		FUNCTION(DEPARTURE, 1)

DEFINE_ENUM(TargetOrDeparture, TARGET_OR_DEPARTURE_VALUES)

// Indique quand la bulle est définie comment interpréter ses valeurs
#define BUBBLE_INTERPRETATION_VALUES(FUNCTION)                          \
		/* Les points de la bulle sont les points */                        \
		/* à distance de norme infinie sur chaque axe*/                     \
		FUNCTION(MOORE, 0)                                                  \
		/* Identique à MOORE mais l'unité des rayons de bulles */           \
		/* est en points (pixels) de grille */                              \
		FUNCTION(MOORE_PX, 1)                                               \
		FUNCTION(EUCLIDEAN, 2)                                              \
		/* Les points de la bulle sont tous les points */                   \
		/* à une distance euclidienne du centre supposée */                 \
		/* en points (pixels) de grille */                                  \
		/* (cas particulier de ELLIPTIC_PX */                               \
		/* où tous les rayons de la bulle sont égaux) */                    \
		FUNCTION(EUCLIDEAN_PX, 3)                                           \
		/* Les points de la bulle sont tous les points */                   \
		/* dans l'ellipsoïde partant du centre */                           \
		/* dont les demi-axes valent les valeurs de BUBBLE_RADIUS */        \
		/* supposée en points (pixels) de grille */                         \
		FUNCTION(ELLIPTIC, 4)                                               \
		FUNCTION(ELLIPTIC_PX, 5)                                            \
		/* Fonction utilisateur isValidNeighbor */                          \
		FUNCTION(CUSTOM, 6)                                                 \
		/* Identique à CUSTOM mais l'unité des rayons de bulles */          \
		/* est en points (pixels) de grille*/                               \
		FUNCTION(CUSTOM_PX, 7)

DEFINE_ENUM(BubbleInterpretation, BUBBLE_INTERPRETATION_VALUES)

inline bool isInPixelsUnits(BubbleInterpretation interp)
    {
    return interp == MOORE_PX || interp == EUCLIDEAN_PX || interp == ELLIPTIC_PX || interp == CUSTOM_PX;
    }

// Ces valeurs sont utilisé dans 
#define STRATEGY_VALUES(FUNCTION)               \
		FUNCTION(FIRST, 0)                          \
		FUNCTION(HEAVY, 1)                          \
		FUNCTION(FIRST_ONLY, 2)                     \
		FUNCTION(SHUFFLE, 3)                        \
		FUNCTION(SORT, 4)                           \
		FUNCTION(RESET_ORDER, 5)                    \
		FUNCTION(CLOSEST, 6)                        \
		FUNCTION(PREFERED, 7)                       \
		FUNCTION(BUBBLE, 8)                         \
		FUNCTION(SMOOTH, 9)                         \
		FUNCTION(TEMPORAL_CONTROL, 10)              \
		FUNCTION(BUBBLE_BORDER, 11)                 \
		FUNCTION(USER_STRATEGY, 255)                \

DEFINE_ENUM(PredefinedStrategyName, STRATEGY_VALUES)

#define TYCHE_DISTRIBUTION_VALUES(FUNCTION)              \
		FUNCTION(UNIFORM, 0)                                 \
		FUNCTION(CONSTANT, 1)                                \
		FUNCTION(CUSTOM_DETERMINED, 2)                       \
		FUNCTION(CUMULATIVE_DISTRIBUTION, 3)                 \
		FUNCTION(CONSTANT_CUMULATIVE_DISTRIBUTION, 4)        \
		FUNCTION(PROBABILITY_DENSITY, 5)                     \
		FUNCTION(CONSTANT_PROBABILITY_DENSITY, 6)            \

DEFINE_ENUM(TycheDistribution, TYCHE_DISTRIBUTION_VALUES)

#endif /* ENUMS_H */
