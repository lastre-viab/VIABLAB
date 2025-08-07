#ifndef EITHER_H
#define EITHER_H

template <typename L, typename R>
class Either {
public:
    Either(const L &left) :
        left(left),
        hasLeft(true) {};
    Either(const R &right) :
        right(right),
        hasLeft(false) {};
    
    template <typename F>    
    const Either &applyLeft(const F &functor) const {
        const Either *res = this;
        if (hasLeft) {
            *res = Either(functor(left));
        }
        return *res;
    }

    template <typename F>    
    Either &applyLeft(const F &functor) {
        Either *res = this;
        if (hasLeft) {
            *res = Either(functor(left));
        }
        return *res;
    }

    bool operator==(const Either &other) const {
        if (isLeft()) {
            return other.isLeft() && left == other.left;
        }
        else {
            return other.isRight() && right == other.right;
        }
    }

    bool operator==(const L &other) const {
        return isLeft() && left == other;
    }

    bool operator==(const R &other) const {
        return isRight() && right == other;
    }

    bool operator!=(const Either &other) const {
        return !this->operator==(other);
    }

    bool operator!=(const L&other) const {
        return !this->operator==(other);
    }

    bool operator!=(const R&other) const {
        return !this->operator==(other);
    }
    
    bool isLeft() const {
        return hasLeft;
    }
    bool isRight() const {
        return !isLeft();
    }
    
    L fromLeft(L defaultValue) const {
        if (isLeft()) {
            return left;
        }
        else {
            return defaultValue;
        }
    }
    R fromRight(R defaultValue) const {
        if (isRight()) {
            return right;
        }
        else {
            return defaultValue;
        }
    }
private:
    union {
        L left;
        R right;
    };
    bool hasLeft;
};

#endif /* EITHER_H */
