/*== SAGITTARIUS =====================================================================
 * Copyright (c) 2012, Jesse K Medley
 * All rights reserved.

 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of The University of Washington nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

//== BEGINNING OF CODE ===============================================================

//== INCLUDES ========================================================================

#include "sbnw/error.h"
#include "network.h"
#include "sbnw/io.h"
#include "sbnw/dist.h"
#include "sbnw/rand_unif.h"
#include "sbnw/allen.h"
#include "sbnw/sign_mag.h"
#include "sbnw/geom.h"

#include <exception>
#include <typeinfo>
#include <cmath>
#include <cstdlib> //rand

static std::string default_comp_id_;
using namespace sbnw;

namespace sbnw {

    // Does elt stand for element? Who knows? Comments would be handy here.
    std::string eltTypeToStr(const NetworkEltType t) {
        switch (t) {
            case NET_ELT_TYPE_SPEC:
                return "Type Species";
            case NET_ELT_TYPE_RXN:
                return "Type Reaction";
            case NET_ELT_TYPE_COMP:
                return "Type Compartment";
            default:
                throw std::invalid_argument("Invalid Type");
        }
    }

    void dumpEltType(std::ostream &os, const NetworkEltType t, uint32_t ind) {
        os << eltTypeToStr(t);
    }

    bool haveDefaultCompartmentId() {
        return !default_comp_id_.empty();
    }

    void setDefaultCompartmentId(const std::string &id) {
        default_comp_id_ = id;
    }

    std::string getDefaultCompartmentId() {
        return default_comp_id_;
    }

    RxnRoleType SBMLRole2GraphfabRole(libsbml::SpeciesReferenceRole_t &role) {
        switch (role) {
            case libsbml::SPECIES_ROLE_SUBSTRATE:
                return RXN_ROLE_SUBSTRATE;
            case libsbml::SPECIES_ROLE_PRODUCT:
                return RXN_ROLE_PRODUCT;
            case libsbml::SPECIES_ROLE_SIDESUBSTRATE:
                return RXN_ROLE_SIDESUBSTRATE;
            case libsbml::SPECIES_ROLE_SIDEPRODUCT:
                return RXN_ROLE_SIDEPRODUCT;
            case libsbml::SPECIES_ROLE_MODIFIER:
                return RXN_ROLE_MODIFIER;
            case libsbml::SPECIES_ROLE_ACTIVATOR:
                return RXN_ROLE_ACTIVATOR;
            case libsbml::SPECIES_ROLE_INHIBITOR:
                return RXN_ROLE_INHIBITOR;
            case libsbml::SPECIES_ROLE_UNDEFINED:
                throw std::invalid_argument("Cannot convert role SPECIES_ROLE_UNDEFINED");
            default:
                throw std::invalid_argument("Unknown role");
                return RXN_ROLE_SUBSTRATE;
        }
    }

    // true for substrates/products
    static bool isRoleActive(RxnRoleType role) {
        switch (role) {
            case RXN_ROLE_SUBSTRATE:
            case RXN_ROLE_PRODUCT:
            case RXN_ROLE_SIDESUBSTRATE:
            case RXN_ROLE_SIDEPRODUCT:
                return true;
            default:
                return false;
        }
    }

    // modifiers match activators and inhibitors
    static bool isGenericModifier(RxnRoleType role) {
        return role == RXN_ROLE_MODIFIER || role == RXN_ROLE_ACTIVATOR || role == RXN_ROLE_INHIBITOR;
    }

    // modifiers match activators and inhibitors
    static bool matchSBML_RoleGenericMod(RxnRoleType u, RxnRoleType v) {
        if (isGenericModifier(u) && isGenericModifier(v))
            return true;
        else
            return u == v;
    }

    //-- Reaction Curves --

    ArrowheadStyle SubCurve::getArrowheadStyle() const {
        return ArrowheadStyleLookup(this);
    }

    ArrowheadStyle PrdCurve::getArrowheadStyle() const {
        return ArrowheadStyleLookup(this);
    }

    ArrowheadStyle ActCurve::getArrowheadStyle() const {
        return ArrowheadStyleLookup(this);
    }

    ArrowheadStyle InhCurve::getArrowheadStyle() const {
        return ArrowheadStyleLookup(this);
    }

    ArrowheadStyle ModCurve::getArrowheadStyle() const {
        return ArrowheadStyleLookup(this);
    }

    //--CLASS NetworkElement--

    void NetworkElement::resetActivity() {
        _v = Point(0., 0.);
    }

    void NetworkElement::doMotion(const double scale) {
        if (_lock)
            return;
        if (_type != NET_ELT_TYPE_COMP)
            throw std::logic_error("bad");
        if (_v.mag2() > 1e-6)
            _p = _p + _v.normed() * scale;
    }

    void NetworkElement::addDelta(const Point &d) {
        _v = _v + d;
    }

    void NetworkElement::capDelta(const double cap) {
        _v = _v.capMag(cap);
    }

    void NetworkElement::capDelta2(const double cap2) {
        _v.capMag2_(cap2);
    }

    void NetworkElement::setCentroid(const Point &p) {
        _p = p;
        _pset = 1;
        recalcExtents();
    }

    void NetworkElement::setGlobalCentroid(const Point &p) {
        _p = itf_ * p;
        _pset = 1;
        recalcExtents();
    }

    Point NetworkElement::getCentroid(COORD_SYSTEM coord) const {
        if (coord == COORD_SYSTEM_LOCAL)
            return _p;
        else if (coord == COORD_SYSTEM_GLOBAL)
            return tf_ * _p;
        else {
            throw std::logic_error("Unknown coord system");
            return _p;
        }
    }

    double NetworkElement::distance(const NetworkElement &e) const {
        NetworkEltShape sp1 = getShape(), sp2 = e.getShape();

        if (sp1 == sp2 && sp2 == ELT_SHAPE_ROUND) {
            double r = euclidean2d(getCentroid(), e.getCentroid()) - radius() - e.radius();
            return std::max(r, 0.);
        } else {
            // works for boxes & mixed (well enough)
            double u = allenDist(getMinX(), getMaxX(), e.getMinX(), e.getMaxX());
            double v = allenDist(getMinY(), getMaxY(), e.getMinY(), e.getMaxY());
            return sqrt(u * u + v * v);
        }
    }

    bool NetworkElement::overlap(const NetworkElement &e) const {
        return (distance(e) == 0.);
    }

    Point NetworkElement::forceVec(const NetworkElement &e) const {
        NetworkEltShape sp1 = getShape(), sp2 = e.getShape();

        if (sp1 == sp2 && sp2 == ELT_SHAPE_ROUND) {
            return (getCentroid() - e.getCentroid()).normed();
        } else {
            if (overlap(e)) {
                //repel via centroids when elements are overlapping
                return (getCentroid() - e.getCentroid()).normed();
            }
            // works for boxes & mixed (well enough)
            double u = -allenOrdered(getMinX(), getMaxX(), e.getMinX(), e.getMaxX());
            double v = -allenOrdered(getMinY(), getMaxY(), e.getMinY(), e.getMaxY());
            return Point(u, v).normed();
        }
    }

    Point NetworkElement::centroidDisplacementFrom(const NetworkElement &e) const {
        return getCentroid() - e.getCentroid();
    }

    void NetworkElement::forceVec_(const NetworkElement &e, Point &p) const {
        NetworkEltShape sp1 = getShape(), sp2 = e.getShape();

        if (sp1 == sp2 && sp2 == ELT_SHAPE_ROUND) {
            p = (getCentroid() - e.getCentroid()).normed();
        } else {
            if (overlap(e)) {
                //repel via centroids when elements are overlapping
                p = (getCentroid() - e.getCentroid()).normed();
                return;
            }
            // works for boxes & mixed (well enough)
            double u = -allenOrdered(getMinX(), getMaxX(), e.getMinX(), e.getMaxX());
            double v = -allenOrdered(getMinY(), getMaxY(), e.getMinY(), e.getMaxY());
            p.x = u;
            p.y = v;
            p.norm_();
        }
    }

    //--CLASS Node--
    void Node::setName(const std::string &name) {
        _name = name;
    }

    const std::string &Node::getName() const {
        return _name;
    }

    const std::string &Node::getId() const {
        return _id;
    }

    void Node::setId(const std::string &id) {
        _id = id;
    }

    const std::string &Node::getGlyph() const {
        return _gly;
    }

    void Node::setGlyph(const std::string &id) {
        _gly = id;
    }

    int Node::alias(network *net) {
        if (!net->containsNode(this))
            throw std::logic_error("No such node in network: network::alias");

        net->clearExcludeFromSubgraphEnum();
        int nsub_before = net->getNumSubgraphs();

        net->clearExcludeFromSubgraphEnum();
        setExcludeFromSubgraphEnum();
        int nsub_after = net->getNumSubgraphs();

        if (nsub_before != nsub_after)
            return 1;

        for (auto i = net->RxnsBegin(); i != net->RxnsEnd(); ++i) {
            Reaction *r = *i;
            int k = 0;

            typedef std::vector<std::pair<Reaction *, Node *> > RxnList;

            RxnList rxnlist;

            for (auto ci = r->CurvesBegin(); ci != r->CurvesEnd(); ++ci, ++k) {
                RxnBezier *c = *ci;

                if (c->ns != this && c->ne != this)
                    continue;

                Node *n = new Node();

                n->setName(getName());

                n->setWidth(getWidth());
                n->setHeight(getHeight());

                {
                    std::stringstream ss;
                    ss << getId() << "_alias" << k;
                    n->setGlyph(ss.str());
                }
                n->setId(getId());
                n->numUses() = 1;
                n->setAlias(true);

                if (Compartment *comp = net->findContainingCompartment(this))
                    comp->addElt(n);
                n->set_i(net->getUniqueIndex());

                n->setCentroid(new2ndPos(c->getCentroidCP(), getCentroid(), 0., -50., false));
                n->setTransform(tf_, false);
                n->setInverseTransform(itf_, false);

                net->addNode(n);

                rxnlist.push_back(std::make_pair(r, n));
            }

            // add to reactions
            for (RxnList::iterator j = rxnlist.begin(); j != rxnlist.end(); ++j)
                j->first->addSpeciesRef(j->second, j->first->getSpeciesRole(this));

            // rebuild curves
            for (RxnList::iterator j = rxnlist.begin(); j != rxnlist.end(); ++j)
                j->first->rebuildCurves();

        }

        try {
            net->removeNode(this);
        } catch (...) {
            throw std::logic_error("Could not remove original node: network::alias");
        }

        return 0;
    }

    bool Node::isCommonInstance(const Node *other) const {
        return getId() == other->getId();
    }

    Point Node::getUpperLeftCorner() const {
        return _p - Point(40, 20);
    }

    Point Node::getLowerRightCorner() const {
        return _p + Point(40, 20);
    }

    void Node::setWidth(double w) {
        Point d(w / 2., getHeight() / 2.);
        _ext.setMin(getCentroid() - d);
        _ext.setMax(getCentroid() + d);
    }

    void Node::setHeight(double h) {
        Point d(getWidth() / 2., h / 2.);
        _ext.setMin(getCentroid() - d);
        _ext.setMax(getCentroid() + d);
    }

    void Node::affectGlobalWidth(double ww) {
        double w = ww / tf_.scaleFactor();
        Point d(w / 2., getHeight() / 2.);
        _ext.setMin(getCentroid() - d);
        _ext.setMax(getCentroid() + d);
    }

    void Node::affectGlobalHeight(double hh) {
        double h = hh / tf_.scaleFactor();
        Point d(getWidth() / 2., h / 2.);
        _ext.setMin(getCentroid() - d);
        _ext.setMax(getCentroid() + d);
    }

    void Node::dump(std::ostream &os, uint32_t ind) {
        indent(os, ind);
        if (isAlias())
            os << "Alias ";
        os << "Node:\n";
        indent(os, ind + 2);
        os << "Name: \"" << _name << "\"\n";
        indent(os, ind + 2);
        os << "ID: \"" << _id << "\"\n";
        if (_comp) {
            indent(os, ind + 2);
            os << "Compartment: " << _comp->getId() << "\n";
        }
        indent(os, ind + 2);
        os << "Degree: " << _deg << "\n";
        indent(os, ind + 2);
        os << "Local degree: " << _ldeg << "\n";
        indent(os, ind + 2);
        os << "Glyph: \"" << _gly << "\"\n";
        indent(os, ind + 2);
        os << "Bounding Box: " << getUpperLeftCorner() << ", " << getLowerRightCorner() << "\n";
    }

    void Node::dumpForces(std::ostream &os, uint32_t ind) const {
        indent(os, ind);
        os << "Node forces: " << _v << "\n";
    }

    std::string rxnRoleToString(RxnRoleType role) {
        switch (role) {
            case RXN_ROLE_SUBSTRATE:
                return "substrate";
            case RXN_ROLE_PRODUCT:
                return "product";
            case RXN_ROLE_SIDESUBSTRATE:
                return "side substrate";
            case RXN_ROLE_SIDEPRODUCT:
                return "side product";
            case RXN_ROLE_MODIFIER:
                return "modifier";
            case RXN_ROLE_ACTIVATOR:
                return "activator";
            case RXN_ROLE_INHIBITOR:
                return "inhibitor";
        }
    }

    std::string CurveTypeToString(RxnCurveType t) {
        switch (t) {
            case RXN_CURVE_SUBSTRATE:
                return "Substrate";
            case RXN_CURVE_PRODUCT:
                return "Product";
            case RXN_CURVE_MODIFIER:
                return "Modifier";
            case RXN_CURVE_ACTIVATOR:
                return "Activator";
            case RXN_CURVE_INHIBITOR:
                return "Inhibitor";
            default:
                return "Unknown";
        }
    }

    // -- CLASS RxnCurveFactory

    RxnBezier *RxnCurveFactory::CreateCurve(RxnRoleType role) {
        switch (role) {
            case RXN_ROLE_SUBSTRATE:
            case RXN_ROLE_SIDESUBSTRATE:
                return new SubCurve();
            case RXN_ROLE_PRODUCT:
            case RXN_ROLE_SIDEPRODUCT:
                return new PrdCurve();
            case RXN_ROLE_MODIFIER:
                return new ModCurve();
            case RXN_ROLE_ACTIVATOR:
                return new ActCurve();
            case RXN_ROLE_INHIBITOR:
                return new InhCurve();
            default:
                throw std::logic_error("Unrecognized species type");
        }
    }

    //--CLASS Reaction--

    void Reaction::hierarchRelease() {
        deleteCurves();
    }

    void Reaction::addSpeciesRef(Node *n, RxnRoleType role) {
        _spec.push_back(std::make_pair(n, role));
        // recompute curves
        _cdirty = true;
        // increase degree
        ++_deg;
        ++_ldeg;
        ++n->_deg;
        ++n->_ldeg;
    }

    void Reaction::removeNode(Node *n) {
        bool rebuild = false;
        repeat:
        for (auto i = _spec.begin(); i != _spec.end(); ++i) {
            Node *x = i->first;
            if (x == n) {
                rebuild = true;
                std::cout << "Rxn: element erased\n";
                --_deg;
                --_ldeg;
                --n->_deg;
                --n->_ldeg;
                _spec.erase(i);
                goto repeat; // in case the species shows up multiple times
            }
        }
        if (rebuild)
            rebuildCurves();
    }

    Node *Reaction::findSpeciesById(const std::string &id) {
        for (auto i = _spec.begin(); i != _spec.end(); ++i) {
            Node *n = i->first;
            if (n->getId() == id)
                return n;
        }
        //not found
        return nullptr;
    }

    bool Reaction::hasSpecies(const Node *n) const {
        for (auto i = NodesBegin(); i != NodesEnd(); ++i) {
            const Node *nn = i->first;
            if (nn == n)
                return true;
        }
        //not found
        return false;
    }

    uint64_t Reaction::degree(const Node *n) {
        unsigned long result = 0;
        for (ConstNodeIt i = NodesBegin(); i != NodesEnd(); ++i) {
            const Node *nn = i->first;
            if (nn == n)
                ++result;
        }
        return result;
    }

    void Reaction::substituteSpeciesById(const std::string &id, Node *spec) {
        for (auto i = _spec.begin(); i != _spec.end(); ++i) {
            Node *n = i->first;
            if (n->getId() == id) {
                --n->_ldeg;
                ++spec->_ldeg;
                i->first = spec;
            }
        }
    }

    void Reaction::substituteSpeciesByIdwRole(const std::string &id, Node *spec, RxnRoleType role) {
        for (auto i = _spec.begin(); i != _spec.end(); ++i) {
            Node *n = i->first;
            if (n->getId() == id && matchSBML_RoleGenericMod(i->second, role)) {
                --n->_ldeg;
                ++spec->_ldeg;
                i->first = spec;
                // SBML inconsistency
                if ((i->second == RXN_ROLE_MODIFIER) && (role == RXN_ROLE_ACTIVATOR || role == RXN_ROLE_INHIBITOR)) {
//                   std::cerr << "Set role for " << spec->getId() << " to " << rxnRoleToString(role) << "\n";
                    i->second = role;
                }
            }
        }
    }

    RxnRoleType Reaction::getSpeciesRole(Node *x) {
        for (auto i = _spec.begin(); i != _spec.end(); ++i) {
            Node *n = i->first;
            if (n == x) {
                return i->second;
            }
        }
        throw std::logic_error("No such node: Reaction::getSpeciesRole");
    }

    void Reaction::substituteSpecies(Node *before, Node *after) {
        for (auto i = _spec.begin(); i != _spec.end(); ++i) {
            Node *n = i->first;
            if (n == before) {
                --n->_ldeg;
                ++after->_ldeg;
                i->first = after;
            }
        }
    }

    Reaction::CurveVec &Reaction::getCurves() {
        curveGuard();
        return _curv;
    }

#define REBUILD_CURVES_DIAG 0

    void Reaction::rebuildCurves() {
        deleteCurves();

# if REBUILD_CURVES_DIAG
        std::cerr << "Rebuild curves\n";
# endif

        for (ConstNodeIt i = NodesBegin(); i != NodesEnd(); ++i) {
            Node *n = i->first;
            RxnRoleType r = i->second;
            // the curve
            RxnBezier *curv = nullptr;
# if REBUILD_CURVES_DIAG
            std::cerr << "  Role: " << rxnRoleToString(r) << "\n";
# endif
            switch (r) {
                case RXN_ROLE_SUBSTRATE:
                case RXN_ROLE_SIDESUBSTRATE:
                    curv = new SubCurve();
                    curv->as = &n->_p;
                    curv->ns = n;
                    curv->owns = 0; //weak ref
                    curv->ae = &_p;
                    curv->owne = 0; //weak ref
                    break;
                case RXN_ROLE_PRODUCT:
                case RXN_ROLE_SIDEPRODUCT:
                    curv = new PrdCurve();
                    curv->as = &_p;
                    curv->owns = 0; //weak ref
                    curv->ae = &n->_p;
                    curv->ne = n;
                    curv->owne = 0; //weak ref
                    break;
                case RXN_ROLE_MODIFIER:
                    curv = new ModCurve();
                    curv->as = &n->_p;
                    curv->ns = n;
                    curv->owns = 0; //weak ref
                    curv->ae = &_p;
                    curv->owne = 0; //weak ref
                    break;
                case RXN_ROLE_ACTIVATOR:
                    curv = new ActCurve();
                    curv->as = &n->_p;
                    curv->ns = n;
                    curv->owns = 0; //weak ref
                    curv->ae = &_p;
                    curv->owne = 0; //weak ref
                    break;
                case RXN_ROLE_INHIBITOR:
                    curv = new InhCurve();
                    curv->as = &n->_p;
                    curv->ns = n;
                    curv->owns = 0; //weak ref
                    curv->ae = &_p;
                    curv->owne = 0; //weak ref
                    break;
                default:
                    std::cerr << "Unrecognized role type\n";
                    throw std::logic_error("Unrecognized role type");
            }
            curv->setTransform(tf_);
            curv->setInverseTransform(itf_);
            _curv.push_back(curv);
        }

        recalcCurveCPs();

        _cdirty = false;

    }

    void Reaction::recalcCurveCPs() {
        uint64_t csub = 0;
        Point ctrlCent(0, 0);
        Point loopPt;
        bool looped = false;

        for (ConstNodeIt i = NodesBegin(); i != NodesEnd(); ++i) {
            Node *n = i->first;
            RxnRoleType r = i->second;

            switch (r) {
                case RXN_ROLE_SUBSTRATE:
                case RXN_ROLE_SIDESUBSTRATE:
                    // control pt stuff
                    ctrlCent += n->getCentroid();
                    csub++;
                    for (ConstNodeIt j = NodesBegin(); j != NodesEnd(); ++j) {
                        Node *nn = j->first;
                        if (nn == n && r != j->second) {
                            looped = true;
                            loopPt = nn->getCentroid();
                        }
                    }
                    break;
                case RXN_ROLE_PRODUCT:
                case RXN_ROLE_SIDEPRODUCT:
                case RXN_ROLE_MODIFIER:
                case RXN_ROLE_ACTIVATOR:
                case RXN_ROLE_INHIBITOR:
                    break;
                default:
                    std::cerr << "Unrecognized role type2\n";
                    throw std::invalid_argument("Unrecognized species type");
            }

        }

        ctrlCent = (ctrlCent + _p) * (1. / (csub + 1));
        double scalar = 20.;

        if (looped) {
            const double d = -scalar;
            ctrlCent = _p + (_p - loopPt);

            ctrlCent = new2ndPos(loopPt, _p, 0., d, false);

            ctrlCent = new2ndPos(_p, ctrlCent, -90., 0., false);
        }

        // Correction applied to uni-uni reactions
        if (NetworkElement::degree() == 2) {
            double d = -(_p - ctrlCent).mag();
            Point p1, p2;

            for (ConstNodeIt i = NodesBegin(); i != NodesEnd(); ++i) {
                Node *n = i->first;
                RxnRoleType r = i->second;

                switch (r) {
                    case RXN_ROLE_SUBSTRATE:
                    case RXN_ROLE_SIDESUBSTRATE:
                        p2 = n->getMin();
                        break;
                    case RXN_ROLE_PRODUCT:
                    case RXN_ROLE_SIDEPRODUCT:
                        p1 = n->getMin();
                        break;
                    case RXN_ROLE_MODIFIER:
                    case RXN_ROLE_ACTIVATOR:
                    case RXN_ROLE_INHIBITOR:
                        break;
                    default:
                        throw std::invalid_argument("Unrecognized species type");
                }
            }

            ctrlCent = _p + (p2 - p1);
            ctrlCent = new2ndPos(ctrlCent, _p, 0, d, false);
        }

        // keep dir, subtract 25 from length
        ctrlCent = new2ndPos(ctrlCent, _p, 0., -scalar, false);

        // control points
        for (auto i = CurvesBegin(); i != CurvesEnd(); ++i) {
            RxnBezier *c = *i;
            RxnCurveType role = c->getRole();

            Box bs(*c->as - Point(scalar * 3 / 2, scalar), *c->as + Point(scalar * 3 / 2, scalar));
            Box be(*c->ae - Point(scalar * 3 / 2, scalar), *c->ae + Point(scalar * 3 / 2, scalar));

            switch (role) {
                case RXN_CURVE_SUBSTRATE:
                    c->s = calcCurveBackup(ctrlCent, *c->as, c->ns ? c->ns->getBoundingBox() : bs, scalar / 2);
                    c->c1 = new2ndPos(_p, c->s, 0., -scalar, false);
                    c->e = *c->ae;
                    c->c2 = ctrlCent;
                    break;
                case RXN_CURVE_PRODUCT:
                    c->s = *c->as;
                    c->c1 = new2ndPos(ctrlCent, _p, 0., 1., true);
                    c->e = calcCurveBackup(c->c1, *c->ae, c->ne ? c->ne->getBoundingBox() : be, scalar / 2);
                    c->c2 = new2ndPos(_p, c->e, 0., -scalar, false);
                    break;
                case RXN_CURVE_ACTIVATOR:
                case RXN_CURVE_INHIBITOR:
                case RXN_CURVE_MODIFIER:
                    c->s = calcCurveBackup(_p, *c->as, c->ns ? c->ns->getBoundingBox() : bs, scalar / 2);
                    c->c1 = new2ndPos(*c->as, _p, 0., -15., false);
                    c->e = c->c1;
                    c->c2 = new2ndPos(*c->as, _p, 0., -20., false);
                    break;
                default:
                    c->s = calcCurveBackup(_p, *c->as, c->ns ? c->ns->getBoundingBox() : bs, scalar / 2);
                    c->c1 = c->s;
                    c->e = _p;
                    c->c2 = _p;
                    break;
            }
        }

        int k_i = 0;
        for (auto i = CurvesBegin(); i != CurvesEnd(); ++i) {
            RxnBezier *c1 = *i;

            auto j = i;
            ++j;
            for (; j != CurvesEnd(); ++j) {
                RxnBezier *c2 = *j;

                if (c1->getNodeUsed() != nullptr && c1->getNodeUsed() == c2->getNodeUsed() &&
                    c1->getRole() == c2->getRole()) {

                    c1->setNodeSideCP(
                            new2ndPos(c1->getNodeUsed()->getCentroid(), c1->getNodeSideCP(), scalar, scalar / 2,
                                      false));
                    c2->setNodeSideCP(
                            new2ndPos(c2->getNodeUsed()->getCentroid(), c2->getNodeSideCP(), -scalar, scalar / 2,
                                      false));

                    c1->setNodeSide(new2ndPos(c1->getNodeSideCP(), c1->getNodeSide(), -scalar, 0., false));
                    c2->setNodeSide(new2ndPos(c2->getNodeSideCP(), c2->getNodeSide(), scalar, 0., false));
                }
            }
        }
    }

    void Reaction::clipCurves(const double padding, const double clip_cutoff) {
        // control points
        for (auto i = CurvesBegin(); i != CurvesEnd(); ++i) {
            RxnBezier *c = *i;
            RxnCurveType role = c->getRole();

            switch (role) {
                case RXN_CURVE_SUBSTRATE:
                case RXN_CURVE_ACTIVATOR:
                case RXN_CURVE_INHIBITOR:
                case RXN_CURVE_MODIFIER:
                    if (c->ns) {
                        Box b = c->ns->getBoundingBox().padded(padding);
                        c->clipReverseToBox(b, clip_cutoff);
                    }
                    break;
                case RXN_CURVE_PRODUCT:
                    if (c->ne) {
                        Box b = c->ne->getBoundingBox().padded(padding);
                        c->clipForwardToBox(b, clip_cutoff);
                    }
                    break;
                default:
                    throw std::invalid_argument("Unrecognized curve type");
                    break;
            }
        }
    }

    void Reaction::recenter() {
//         std::cerr << "RECENTER\n";
        uint32_t count = 0;
        _p = Point(0., 0.);
        for (ConstNodeIt i = NodesBegin(); i != NodesEnd(); ++i) {
            Node *n = i->first;
            _p = _p + n->getCentroid();
            ++count;
        }
        // normalize
        _p = _p * (1. / count);
        rebuildCurves();
    }

    void Reaction::recompCentroid() {
//         std::cerr << "RECOMP CENTROID\n";
        if (isCentroidSet())
            return;
        doCentroidCalc();
    }

    void Reaction::forceRecalcCentroid() {
//         std::cerr << "RECALC CENTROID\n";
        doCentroidCalc();
        _pset = 1;
    }

    void Reaction::doCentroidCalc() {
        uint32_t count = 0;
        _p = Point(0., 0.);
        for (ConstNodeIt i = NodesBegin(); i != NodesEnd(); ++i) {
            // detect duplicates
            for (ConstNodeIt j = NodesBegin(); j != i; ++j)
                if (i->first == j->first)
                    goto doCentroidCalc_skip;

            {
                Node *n = i->first;
                _p = _p + n->getCentroid();
                ++count;
            }

            doCentroidCalc_skip:;
        }
        // normalize
        _p = _p * (1. / count);
    }

    void Reaction::deleteCurves() {
        for (auto i = _curv.begin(); i != _curv.end(); ++i) {
            delete *i;
        }
        _curv.clear();
    }

    void Reaction::dump(std::ostream &os, uint32_t ind) {
        indent(os, ind);
        os << "Reaction:\n";
        indent(os, ind + 2);
        os << "ID: \"" << _id << "\"\n";
        indent(os, ind + 2);
        os << "Degree: " << _deg << "\n";
        indent(os, ind + 2);
        os << "Local degree: " << _ldeg << "\n";
        indent(os, ind + 2);
        os << "Species: \n";
        for (ConstNodeIt i = _spec.begin(); i != _spec.end(); ++i) {
            indent(os, ind + 4);
            os << i->first->getId() << "(" << i->first->getGlyph() << "), role: " << rxnRoleToString(i->second) << "\n";
        }
    }

    void Reaction::dumpForces(std::ostream &os, uint32_t ind) const {
        indent(os, ind);
        os << "Reaction forces: " << _v << "\n";
    }

    //--CLASS Compartment--

    void Compartment::addElt(NetworkElement *e) {
        _elt.push_back(e);
    }

    bool Compartment::containsElt(const NetworkElement *e) const {
        for (auto i = EltsBegin(); i != EltsEnd(); ++i) {
            const NetworkElement *x = *i;
            if (x == e)
                return true;
        }
        return false;
    }

    void Compartment::removeElt(NetworkElement *e) {
        for (auto i = EltsBegin(); i != EltsEnd(); ++i) {
            if (*i == e) {
                printf("Element erased\n");
                _elt.erase(i);
                return;
            }
        }
    }

    void Compartment::setRestExtents(const Box &ext) {
        _ext = ext;
        _ra = _ext.area();
    }

    void Compartment::resizeEnclose(double padding) {
        double minx = 0, miny = 0, maxx = 0, maxy = 0;
        auto i = EltsBegin();
        if (i != EltsEnd()) {
            NetworkElement *e = *i;
            minx = e->getMinX();
            miny = e->getMinY();
            maxx = e->getMaxX();
            maxy = e->getMaxY();
            ++i;
        }
        for (; i != EltsEnd(); ++i) {
            NetworkElement *e = *i;
            minx = std::min(minx, e->getMinX());
            maxx = std::max(maxx, e->getMaxX());
            miny = std::min(miny, e->getMinY());
            maxy = std::max(maxy, e->getMaxY());
        }
        _ext = Box(Point(minx, miny), Point(maxx, maxy));
        _ext = _ext.padded(padding);
        _ra = _ext.area();
    }

    void Compartment::autoSize() {
        uint64_t count = _elt.size();
        double dim = 350 * sqrt((double) count);
        // avoid singularities in layout algo
        Point shake((rand() % 1000) / 100., (rand() % 1000) / 100.);
        _ext = Box(Point(0, 0) + shake, Point(dim, dim) + shake);
        //_ext = Box(Point(0, 0), Point(dim, dim));
        _ra = _ext.area();
    }

    void Compartment::resetActivity() {
        _v = Point(0, 0);
        // now calculate stress due to being stretched beyond rest area
        // this stress always acts to shrink the comp
        double w = _ext.width(), h = _ext.height();
        double d2 = _ext.area() - _ra;
        // strain (liberally speaking), evenly distributed along all axes
        double strain = sign(d2) * sqrt(mag(d2) / _ra);
        _fx1 = _res * _E * strain * w;
        _fy1 = _res * _E * strain * h;
        _fx2 = -_res * _E * strain * w;
        _fy2 = -_res * _E * strain * h;
    }

    void Compartment::applyBoundaryForce(const double fx1, const double fy1, const double fx2, const double fy2) {
        _fx1 += fx1;
        _fy1 += fy1;
        _fx2 += fx2;
        _fy2 += fy2;
    }

    void Compartment::doInternalForce(NetworkElement *e, const double f, const double t) {
        double x1 = _ext.getMin().x, y1 = _ext.getMin().y, x2 = _ext.getMax().x, y2 = _ext.getMax().y;
        double invt = 1. / t;

        double eminx = e->getMinX();
        double eminy = e->getMinY();
        double emaxx = e->getMaxX();
        double emaxy = e->getMaxY();

        // compute forces
        double fx1 = f * exp((x1 - eminx) * invt);
        double fx2 = -f * exp((emaxx - x2) * invt);
        double fy1 = f * exp((y1 - eminy) * invt);
        double fy2 = -f * exp((emaxy - y2) * invt);

        // do forces on element
        e->addDelta(Point(fx1 + fx2, fy1 + fy2));

        // do forces on container
        applyBoundaryForce(-fx1, -fx2, -fy1, -fy2);
        addDelta(-Point(fx1 + fx2, fy1 + fy2));
    }

    void Compartment::doInternalForceAll(const double f, const double t) {
        for (auto i = EltsBegin(); i != EltsEnd(); ++i) {
            NetworkElement *e = *i;
            doInternalForce(e, f, t);
        }
    }

    void Compartment::doMotion(const double scale_) {
        if (_lock)
            return;
        const double scale = 0.2 * scale_;
        double w = _ext.width(), h = _ext.height();
        // adjust the extents based on Hooke's law of elasticity
        // forces -> stress -> strain -> displacement
        _ext.setMin(_ext.getMin() + (scale / _E) * Point(_fx1 * w / h, _fy1 * h / w) + scale * _v);
        _ext.setMax(_ext.getMax() + (scale / _E) * Point(_fx2 * w / h, _fy2 * h / w) + scale * _v);
        if (_ext.width() < 10.)
            _ext.setWidth(10.);
        if (_ext.height() < 10.)
            _ext.setHeight(10.);
        //recalc centroid?
    }

    void Compartment::capDelta2(const double cap2) {
        _v.capMag2_(cap2);
        const double cap = sqrt(cap2);
        if (mag(_fx1) > cap)
            _fx1 = sign(_fx1) * cap;
        if (mag(_fy1) > cap)
            _fy1 = sign(_fy1) * cap;
        if (mag(_fx2) > cap)
            _fx2 = sign(_fx2) * cap;
        if (mag(_fy2) > cap)
            _fy2 = sign(_fy2) * cap;
    }

    bool Compartment::contains(const NetworkElement *e) const {
        for (auto i = EltsBegin(); i != EltsEnd(); ++i) {
            if (*i == e)
                return true;
        }
        return false;
    }

    void Compartment::dump(std::ostream &os, uint32_t ind) {
        indent(os, ind);
        os << "Compartment:\n";
        indent(os, ind + 2);
        os << "ID: \"" << _id << "\"\n";
        indent(os, ind + 2);
        os << "Glyph: \"" << _gly << "\"\n";
        indent(os, ind + 2);
        os << "Extents: " << _ext << "\n";
    }

    void Compartment::dumpForces(std::ostream &os, uint32_t ind) const {
        indent(os, ind);
        os << "Compartment forces: " << _fx1 << ", " << _fy1 << ", " << _fx2 << ", " << _fy2 << "), Centroid forces: "
           << _v << "\n";
    }

    //--CLASS network--

    void network::hierarchRelease() {
        // FIXME: replace with hierarch free
        for (auto i = _nodes.begin(); i != _nodes.end(); ++i) {
            // (*i)->hierarchRelease();
            delete *i;
        }
        for (auto &i : _rxn) {
            i->hierarchRelease();
            delete i;
        }
        for (auto &i : _comp) {
            // (*i)->hierarchRelease();
            delete i;
        }
    }

    void network::addNode(Node *n) {
        if (!n)
            throw std::logic_error("No node to add");
        _nodes.push_back(n);
        addElt(n);
    }

    void network::removeReactionsForNode(Node *n) {
        for (auto i = _rxn.begin(); i != _rxn.end(); ++i) {
            (*i)->removeNode(n);
        }
    }

    void network::removeNode(Node *n) {
        if (!n)
            throw std::logic_error("No node to remove");
        // remove from element container
        removeElt(n);
        // remove from compartments
        for (auto i = CompsBegin(); i != CompsEnd(); ++i) {
            Compartment *c = *i;
            c->removeElt(n);
        }
        removeReactionsForNode(n);
        for (auto i = _nodes.begin(); i != _nodes.end(); ++i) {
            Node *x = *i;
            if (x == n) {
                _nodes.erase(i);
                std::cout << "Removed node " << n << "\n";
                return;
            }
        }
        throw std::logic_error("No such node: network::removeNode");
    }

    void network::connectNode(Node *n, Reaction *r, RxnRoleType role) {
        r->addSpeciesRef(n, role);
        r->rebuildCurves();
    }

    bool network::isNodeConnected(Node *n, Reaction *r) const {
        if (!n)
            throw std::logic_error("No node");
        if (!r)
            throw std::logic_error("No reaction");
        return r->hasSpecies(n);
    }

    Node *network::findNodeById(const std::string &id) {
        for (auto i = _nodes.begin(); i != _nodes.end(); ++i) {
            Node *n = *i;
            if (n->getId() == id)
                return n;
        }
        //not found
        return nullptr;
    }

    const Node *network::findNodeById(const std::string &id) const {
        for (auto i = _nodes.begin(); i != _nodes.end(); ++i) {
            const Node *n = *i;
            if (n->getId() == id)
                return n;
        }
        //not found
        return nullptr;
    }

    std::string network::getUniqueId() const {
        std::size_t k = 0;
        std::string id;
        const Node *n = nullptr;

        do {
            ++k;
            std::stringstream ss;
            ss << "Node_" << k;
            id = ss.str();
            std::cout << "Trying " << id << "\n";
        } while (findNodeById(id));

        std::cout << "Unique ID: " << id << "\n";

        return id;
    }

    std::string network::getUniqueGlyphId(const Node &src) const {
        static std::size_t k = 0;
        const Node *n = nullptr;

        ++k;
        std::stringstream ss;
        ss << src.getGlyph() << "_" << k;

        return ss.str();
    }

    std::size_t network::getUniqueIndex() const {
//         std::cout << "getUniqueIndex started\n";
        std::size_t k = 0;

        repeat:
        for (auto i = _nodes.begin(); i != _nodes.end(); ++i) {
            const Node *node = *i;
            if (node->get_i() == k) {
                ++k;
                goto repeat;
            }
        }

        return k;
    }

    Node *network::findNodeByGlyph(const std::string &gly) {
        for (auto i = _nodes.begin(); i != _nodes.end(); ++i) {
            Node *n = *i;
            if (n->getGlyph() == gly)
                return n;
        }
        //not found
        return nullptr;
    }

    Node *network::getUniqueNodeAt(const size_t n) {
        size_t k = 0, a = 1;
        for (auto i = _nodes.begin(); i != _nodes.end(); ++i) {
            Node *x = *i;
            if (k == n)
                return x;
            if (!x->isAlias()) {
                ++k;
                a = 1;
            } else {
                k += a;
                a = 0;
            }
        }
        {
            std::stringstream ss;
            ss << "No unique node with given index " << n << " where number of unique nodes is " << getNumUniqueNodes();
        }
    }

    size_t network::getNumInstances(const Node *u) {
        size_t k = 0;
        for (auto i = _nodes.begin(); i != _nodes.end(); ++i) {
            Node *v = *i;
            if (u->isCommonInstance(v))
                ++k;
        }
        return k;
    }

    Node *network::getInstance(const Node *u, const size_t n) {
        size_t k = 0;
        for (auto i = _nodes.begin(); i != _nodes.end(); ++i) {
            Node *v = *i;
            if (u->isCommonInstance(v)) {
                if (k == n)
                    return v;
                else
                    ++k;
            }
        }
    }

    bool network::containsNode(const Node *n) const {
        for (auto i = _nodes.begin(); i != _nodes.end(); ++i) {
            const Node *x = *i;
            if (x == n)
                return true;
        }
        return false;
    }

    bool network::containsReaction(const Reaction *r) const {
        for (auto i = RxnsBegin(); i != RxnsEnd(); ++i) {
            const Reaction *x = *i;
            if (x == r)
                return true;
        }
        return false;
    }

    network::AttachedRxnList network::getConnectedReactions(const Node *n) {
        AttachedRxnList result;
        for (ConstRxnIt i = RxnsBegin(); i != RxnsEnd(); ++i) {
            Reaction *x = *i;
            if (x->hasSpecies(n))
                result.push_back(x);
        }
        return result;
    }

    network::AttachedCurveList network::getAttachedCurves(const Node *n) {
        AttachedRxnList rxns = getConnectedReactions(n);
        AttachedCurveList result;
        for (auto i = rxns.begin(); i != rxns.end(); ++i) {
            Reaction *r = *i;
            for (auto j = r->CurvesBegin(); j != r->CurvesEnd(); ++j) {
                RxnBezier *c = *j;
                if (c->includes(n))
                    result.push_back(c);
            }
        }
        return result;
    }

    int network::getNumSubgraphs() {
        enumerateSubgraphs();
        return nsub_;
    }

    void network::enumerateSubgraphs() {
        nsub_ = 0;
        loop:
        for (NodeVec::const_iterator i = _nodes.begin(); i != _nodes.end(); ++i) {
            Node *x = *i;
            if (!x->isSetSubgraphIndex()) {
                propagateSubgraphIndex(x, nsub_++);
                goto loop;
            }
        }
    }

    void network::propagateSubgraphIndex(Node *x, int isub) {
        x->setSubgraphIndex(isub);
        for (auto i = _rxn.begin(); i != _rxn.end(); ++i) {
            Reaction *r = *i;
            if (r->hasSpecies(x)) {
                for (auto j = r->NodesBegin(); j != r->NodesEnd(); ++i) {
                    if (!j->first->isSetSubgraphIndex())
                        propagateSubgraphIndex(j->first, isub);
                }
            }
        }
    }

    void network::clearSubgraphInfo() {
        for (NodeVec::const_iterator i = _nodes.begin(); i != _nodes.end(); ++i) {
            Node *x = *i;
            x->clearSubgraphIndex();
        }
    }

    void network::clearExcludeFromSubgraphEnum() {
        for (NodeVec::const_iterator i = _nodes.begin(); i != _nodes.end(); ++i) {
            Node *x = *i;
            x->clearExcludeFromSubgraphEnum();
        }
    }

    Reaction *network::findReactionById(const std::string &id) {
        for (auto i = _rxn.begin(); i != _rxn.end(); ++i) {
            Reaction *r = *i;
            if (r->getId() == id)
                return r;
        }
        //not found
        return nullptr;
    }

    Compartment *network::findCompById(const std::string &id) {
        for (auto i = _comp.begin(); i != _comp.end(); ++i) {
            Compartment *c = *i;
            if (c->getId() == id)
                return c;
        }
        //not found
        return nullptr;
    }

    Compartment *network::findCompByGlyph(const std::string &gly) {
        for (auto i = _comp.begin(); i != _comp.end(); ++i) {
            Compartment *c = *i;
            if (c->getGlyph() == gly)
                return c;
        }
        //not found
        return nullptr;
    }

    void network::resetUsageInfo() {
        for (auto i = NodesBegin(); i != NodesEnd(); ++i) {
            Node *n = *i;
            n->numUses() = 0;
        }
    }

    void network::addReaction(Reaction *rxn) {
        addElt(rxn);
    }

    void network::removeReaction(Reaction *r) {
        if (!r)
            throw std::logic_error("No reaction to remove");
        // remove from element container
        removeElt(r);
        for (auto i = _rxn.begin(); i != _rxn.end(); ++i) {
            Reaction *x = *i;
            if (x == r) {
                _rxn.erase(i);
                std::cout << "Removed reaction " << r << "\n";
                return;
            }
        }
        throw std::logic_error("No such reaction: network::removeReaction");
    }

    void network::elideEmptyComps() {
        // replace in elt vec
        EltVec w;
        for (auto i = EltsBegin(); i != EltsEnd(); ++i) {
            NetworkElement *e = *i;
            if (e->getType() == NET_ELT_TYPE_COMP) {
                auto *c = (sbnw::Compartment *) e;
                if (!c->empty())
                    w.push_back(c);
            } else
                w.push_back(e);
        }
        _elt.swap(w);

        // replace in comp vec & delete empty ones
        CompVec v;
        for (auto i = CompsBegin(); i != CompsEnd(); ++i) {
            Compartment *c = *i;
            if (!c->empty())
                v.push_back(c);
            else
                delete c;
        }
        _comp.swap(v);
    }

    Compartment *network::findContainingCompartment(const NetworkElement *e) {
        for (auto i = CompsBegin(); i != CompsEnd(); ++i) {
            Compartment *c = *i;
            if (c->containsElt(e))
                return c;
        }
        return nullptr;
    }

    uint64_t network::getNumUniqueNodes() const {
        uint64_t k = 0, a = 1;
        for (auto i = _nodes.begin(); i != _nodes.end(); ++i) {
            const Node *x = *i;
            if (!x->isAlias()) {
                ++k;
                a = 1;
            } else {
                k += a;
                a = 0;
            }
        }
        return k;
    }

    Box network::getBoundingBox() const {
        Box b;
        {
            auto i = EltsBegin();
            if (i == EltsEnd())
                return b;
            NetworkElement *e = *i;
            b = e->getBoundingBox();
        }
//         std::cerr << "network initial bounding box: " << b << "\n";
        for (auto i = EltsBegin() + 1; i != EltsEnd(); ++i) {
            NetworkElement *e = *i;
            b.expandx(e->getBoundingBox());
//             std::cerr << "  Expand by: " << e->getBoundingBox() << "\n";
//             std::cerr << "  new bounding box: " << b << "\n";
        }
//         std::cerr << "network bounding box: " << b << "\n";
        return b;
    }

    void network::fitToWindow(const Box &w) {
        sbnw::Affine2d tf = sbnw::Affine2d::FitToWindow(getBoundingBox(), w);
//         std::cerr << "Applying tf:\n" << tf;
        setTransform(tf);
        setInverseTransform(tf.inv());
    }

    void network::applyTransform(const Affine2d &t) {
        for (auto i = EltsBegin(); i != EltsEnd(); ++i) {
            NetworkElement *e = *i;
            e->applyTransform(t);
        }
    }

    void network::setTransform(const Affine2d &t, bool recurse) {
        for (auto i = EltsBegin(); i != EltsEnd(); ++i) {
            NetworkElement *e = *i;
            e->setTransform(t, recurse);
        }
    }

    void network::setInverseTransform(const Affine2d &it, bool recurse) {
        for (auto i = EltsBegin(); i != EltsEnd(); ++i) {
            NetworkElement *e = *i;
            e->setInverseTransform(it, recurse);
        }
    }

    void network::applyDisplacement(const Point &d) {
        for (auto i = EltsBegin(); i != EltsEnd(); ++i) {
            NetworkElement *e = *i;
            e->applyDisplacement(d);
        }
    }

    void network::resetActivity() {
        for (auto i = EltsBegin(); i != EltsEnd(); ++i) {
            NetworkElement *e = *i;
            e->resetActivity();
        }
    }

    void network::updatePositions(const double scale) {
        for (auto i = EltsBegin(); i != EltsEnd(); ++i) {
            NetworkElement *e = *i;
            e->doMotion(scale);
        }
    }

    void network::resizeCompsEnclose(double padding) {
        for (auto i = CompsBegin(); i != CompsEnd(); ++i) {
            Compartment *c = *i;
            c->resizeEnclose(padding);
        }
    }

    void network::autosizeComps() {
        for (auto i = CompsBegin(); i != CompsEnd(); ++i) {
            Compartment *c = *i;
            c->autoSize();
        }
    }

    void network::updateExtents() {
        for (auto i = EltsBegin(); i != EltsEnd(); ++i) {
            NetworkElement *e = *i;
            e->recalcExtents();
        }
    }

//    void capDelta2(const double cap) {
//        for (auto i = EltsBegin(); i != EltsEnd(); ++i) {
//            NetworkElement *e = *i;
//            e->capDelta2(cap * cap);
//        }
//    }

    Point network::pmean() const {
        Point m(0., 0.);
        uint64_t c = 0;
        for (auto i = EltsBegin(); i != EltsEnd(); ++i) {
            NetworkElement *e = *i;
            m = m + e->getCentroid();
            ++c;
        }
        m = m * (1. / c);
        return m;
    }

    Point network::center() const {
        return getExtents().getCenter();
    }

    Box network::getExtents() const {
        if (EltsBegin() == EltsEnd()) return {};
        Box m((*EltsBegin())->getExtents());

        for (auto i = EltsBegin(); i != EltsEnd(); ++i) {
            NetworkElement *e = *i;
            if (e->getMinX() < m.getMinX())
                m.setMinX(e->getMinX());
            if (e->getMinY() < m.getMinY())
                m.setMinY(e->getMinY());
            if (e->getMaxX() > m.getMaxX())
                m.setMaxX(e->getMaxX());
            if (e->getMaxY() > m.getMaxY())
                m.setMaxY(e->getMaxY());
        }

        return m;
    }

    void network::recenter(const Point &p) {
        Point m(pmean());
        Point d = p - m;
        for (auto i = EltsBegin(); i != EltsEnd(); ++i) {
            NetworkElement *e = *i;
            e->setCentroid(e->getCentroid() + d);
        }
    }

    Point network::pvariance() const {
        Point m(pmean());
        Point d(0., 0.);
        uint64_t c = 0;
        for (auto i = EltsBegin(); i != EltsEnd(); ++i) {
            NetworkElement *e = *i;
            d = d + Point(e->getCentroid() - m).squareTerms();
            ++c;
        }
        d = d.sqrtTerms() * (1. / c);
        return d;
    }

    void network::randomizePositions(const Box &b) {
        for (auto i = _nodes.begin(); i != _nodes.end(); ++i) {
            Node *n = *i;
            if (n->isLocked())
                break;
            n->setCentroid(rand_range(b.getMin().x, b.getMax().x),
                           rand_range(b.getMin().y, b.getMax().y));
        }
        for (auto i = _rxn.begin(); i != _rxn.end(); ++i) {
            Reaction *r = *i;
            if (r->isLocked())
                break;
            r->setCentroid(Point(rand_range(b.getMin().x, b.getMax().x),
                                 rand_range(b.getMin().y, b.getMax().y)));
        }
        for (auto i = CompsBegin(); i != CompsEnd(); ++i) {
            sbnw::Compartment *c = *i;
            if (c->isLocked())
                break;
            double d = sqrt(c->restArea());
            Point p(rand_range(b.getMin().x, b.getMax().x),
                    rand_range(b.getMin().y, b.getMax().y));
            Point dim(d, d);
            c->setExtents(Box(p - dim, p + dim));
        }
        recalcCurveCPs();
        //dump(std::cout, 0);
    }

    void network::rebuildCurves() {
        for (auto i = RxnsBegin(); i != RxnsEnd(); ++i) {
            Reaction *r = *i;
            r->rebuildCurves();
        }
        clipCurves();
    }

    void network::recalcCurveCPs() {
        for (auto i = RxnsBegin(); i != RxnsEnd(); ++i) {
            Reaction *r = *i;
            r->recalcCurveCPs();
        }
    }

    void network::clipCurves(const double padding, const double clip_cutoff) {
        for (auto i = RxnsBegin(); i != RxnsEnd(); ++i) {
            Reaction *r = *i;
            r->clipCurves(padding, clip_cutoff);
        }
    }

    void network::recenterJunctions() {
//         std::cerr << "Recenter junctions\n";
        for (auto i = RxnsBegin(); i != RxnsEnd(); ++i) {
            Reaction *r = *i;
            r->recenter();
        }
    }

    // IO/Diagnostics:

    void network::dump(std::ostream &os, uint32_t ind) {
        indent(os, ind);
        os << "network:\n";
        for (ConstEltIt i = EltsBegin(); i != EltsEnd(); ++i) {
            NetworkElement *e = *i;
            e->dump(os, ind + 2);
        }
    }

    void network::dumpEltForces(std::ostream &os, uint32_t ind) const {
        for (auto i = EltsBegin(); i != EltsEnd(); ++i) {
            NetworkElement *e = *i;
            e->dumpForces(os, ind + 2);
        }
    }

    //--GLOBAL--

    network *networkFromLayout(const libsbml::Layout &lay, const libsbml::Model &mod) {
        network *net = networkFromModel(mod);

        // used to compute aliases
        net->resetUsageInfo();

        //add additional information from layout
        // for compartments
        for (unsigned int i = 0; i < lay.getNumCompartmentGlyphs(); ++i) {
            const libsbml::CompartmentGlyph *cg = lay.getCompartmentGlyph(i);

            sbnw::Compartment *c = net->findCompById(cg->getCompartmentId());

            c->setGlyph(cg->getId());

            const libsbml::BoundingBox *bbox = cg->getBoundingBox();

            c->setRestExtents(
                    Box(Point(bbox->x(), bbox->y()), Point(bbox->x() + bbox->width(), bbox->y() + bbox->height())));
        }

        // place elements inside parent compartments
        for (auto z = net->NodesBegin(); z != net->NodesEnd(); ++z) {
            Node *n = *z;
            Compartment *c = net->findContainingCompartment(n);
            if (c)
                n->setCentroid(c->getCentroid());
        }
        for (auto z = net->RxnsBegin(); z != net->RxnsEnd(); ++z) {
            Reaction *r = *z;
            Compartment *c = net->findContainingCompartment(r);
            if (c)
                r->setCentroid(c->getCentroid());
        }

        // for nodes
        for (int i = 0; i < lay.getNumSpeciesGlyphs(); ++i) {
            const libsbml::SpeciesGlyph *sg = lay.getSpeciesGlyph(i);

            Node *n = net->findNodeById(sg->getSpeciesId());
            if (!n)
                throw std::logic_error("No such node exists");

//             std::cerr << "Species glyph: " << sg->getSpeciesId() << "(" << sg->getId() << ")\n";

            //increment usage counter (used to find aliases)
            if (n->numUses() == 0) {
                n->numUses()++;
                n->setGlyph(sg->getId());
            } else {
                //create an alias node
                n->setAlias(true);
                n = new Node(*n);
                n->setGlyph(sg->getId());
                n->_ldeg = 0;
                //add alias node to the network
                net->addNode(n);
            }

            const libsbml::BoundingBox *bb = sg->getBoundingBox();

            n->setCentroid(Point(bb->x() + bb->width() / 2., bb->y() + bb->height() / 2.));
            n->setWidth(bb->width());
            n->setHeight(bb->height());
        }

        // for reactions
        for (int i = 0; i < lay.getNumReactionGlyphs(); ++i) {
            const libsbml::ReactionGlyph *rg = lay.getReactionGlyph(i);
//             std::cerr << "Read ReactionGlyph " << rg->getId() << "\n";

            Reaction *r = net->findReactionById(rg->getReactionId());
            if (!r)
                throw std::logic_error("No such reaction");

            for (int i_spc = 0; i_spc < rg->getNumSpeciesReferenceGlyphs(); ++i_spc) {
                const libsbml::SpeciesReferenceGlyph *srg = rg->getSpeciesReferenceGlyph(i_spc);

                //get the alias
                Node *alias = net->findNodeByGlyph(srg->getSpeciesGlyphId());
                if (!alias)
                    throw std::logic_error("Unable to find alias node");

                libsbml::SpeciesReferenceRole_t libsbml_role = srg->getRole();
                RxnRoleType role = SBMLRole2GraphfabRole(libsbml_role);
                r->substituteSpeciesByIdwRole(srg->getSpeciesReferenceId(), alias, role);

                // delete preexisting curves
                r->deleteCurves();

                libsbml::Curve const *curve = rg->getCurve();
                libsbml::BoundingBox const *sbml_bb = rg->getBoundingBox();


                // first try bounding box (the proper method, which none of the models use)
                if (sbml_bb &&
                    !(sbml_bb->getPosition()->x() == 0 && sbml_bb->getPosition()->y() == 0 &&
                      sbml_bb->getDimensions() && sbml_bb->getDimensions()->getWidth() == 0 &&
                      sbml_bb->getDimensions()->getHeight() == 0)) {

                    double x_offset = 0, y_offset = 0;
                    if (sbml_bb->getDimensions())
                        x_offset = sbml_bb->getDimensions()->getWidth() * 0.5, y_offset =
                                                                                       sbml_bb->getDimensions()->getHeight() *
                                                                                       0.5;
                    r->setCentroid(sbml_bb->getPosition()->x(), sbml_bb->getPosition()->y());
                } else if (curve && curve->getNumCurveSegments() > 0) {
                    // next try using preexisting centroid coords via reaction curve
//                 std::cerr << "manual centroid\n";
                    r->setCentroid(curve->getCurveSegment(0)->getEnd()->x(), curve->getCurveSegment(0)->getEnd()->y());

                    for (unsigned int j = 0; j < rg->getNumSpeciesReferenceGlyphs(); ++j) {
                        libsbml::SpeciesReferenceRole_t r2 = rg->getSpeciesReferenceGlyph(j)->getRole();
                        RxnRoleType r3 = SBMLRole2GraphfabRole(r2);
                        RxnBezier *c = r->addCurve(r3);
//                     std::cerr << "  Created curve with role " << rxnRoleToString(SBMLRole2GraphfabRole(rg->getSpeciesReferenceGlyph(j)->getRole())) << "\n";
                        Node *target = net->findNodeByGlyph(rg->getSpeciesReferenceGlyph(j)->getSpeciesGlyphId());

                        if (c->getRole() == RXN_CURVE_PRODUCT) {
                            c->as = &r->_p;
                            c->ae = &target->_p;
                            c->ne = target;
                        } else {
                            c->as = &target->_p;
                            c->ns = target;
                            c->ae = &r->_p;
                        }
                        c->owne = 0;
                        c->owns = 0;
                    }

                    r->recalcCurveCPs();

                    r->clearDirtyFlag();

                    // Try to fill in the CP data from layout info
                    for (unsigned int j = 0; j < rg->getNumSpeciesReferenceGlyphs(); ++j) {
                        libsbml::SpeciesReferenceGlyph const *srg = rg->getSpeciesReferenceGlyph(j);
                        RxnBezier *c = r->getCurve(j);
                        // use the first curve segment
                        libsbml::Curve const *sr_curve = srg->getCurve();
                        libsbml::LineSegment const *sr_line = sr_curve->getCurveSegment(0);
                        auto const *sr_bez = dynamic_cast< libsbml::CubicBezier const * >(sr_line);
                        if (sr_bez) {
                            //std::cerr << "sr_bez\n";
                            c->s.x = sr_bez->getStart()->x();
                            c->s.y = sr_bez->getStart()->y();
                            c->e.x = sr_bez->getEnd()->x();
                            c->e.y = sr_bez->getEnd()->y();

                            c->c1.x = sr_bez->getBasePoint1()->x();
                            c->c1.y = sr_bez->getBasePoint1()->y();
                            c->c2.x = sr_bez->getBasePoint2()->x();
                            c->c2.y = sr_bez->getBasePoint2()->y();
                        } else if (sr_line) {
                            //std::cerr << "no sr_bez\n";
                            c->s.x = sr_line->getStart()->x();
                            c->s.y = sr_line->getStart()->y();
                            c->e.x = sr_line->getEnd()->x();
                            c->e.y = sr_line->getEnd()->y();

                            c->c1.x = sr_line->getStart()->x();
                            c->c1.y = sr_line->getStart()->y();
                            c->c2.x = sr_line->getEnd()->x();
                            c->c2.y = sr_line->getEnd()->y();
                            //  CPs should be separated from endpoints for endcap orientation
                            Point ctmp = c->c1;
                            c->c1 = 0.9 * c->c1 + 0.1 * c->c2;
                            c->c2 = 0.9 * c->c2 + 0.1 * ctmp;
                        }
                    }
                } else {
                    // poor results

                    // if all else fails average node coords
                    r->forceRecalcCentroid();

                    skip_recalc_centroid:;
                }
            }

            net->setLayoutSpecified(true);

            return net;
        }
    }

    network *networkFromModel(const libsbml::Model &mod) {
        auto *net = new network();

        if (mod.isSetId())
            net->setId(mod.getId());

        // add compartments
        for (int i = 0; i < mod.getNumCompartments(); ++i) {
            const libsbml::Compartment *comp = mod.getCompartment(i);

            // elide "default" compartments based on SBO
            if (comp->isSetSBOTerm() && comp->getSBOTerm() == 410) {
                continue;
            }

            // assume a compartment with the id "default" or "compartment" represents
            // a default, non-visual compartment, so discard it from the model
            if (comp->getId() != "default" && comp->getId() != "compartment" &&
                comp->getId() != "sbnw_default_compartment" &&
                (!haveDefaultCompartmentId() || getDefaultCompartmentId() != comp->getId())) {
                auto *c = new Compartment();

                // set id
                c->setId(comp->getId());

                // add to network
                net->addCompartment(c);
            }
        }

        // add nodes
        //printf("# floating = %d\n", floating);
        for (int i = 0; i < mod.getNumSpecies(); ++i) {
            Node *n = new Node();

            const libsbml::Species *s = mod.getSpecies(i);

            if (!s)
                throw std::logic_error("Failed to get species");

            n->setName(s->getName());
            n->setId(s->getId());

            //no alias info
            n->numUses() = 1;
            n->setAlias(false);

            // associate compartment (if one exists)
            sbnw::Compartment *c = net->findCompById(s->getCompartment());
            if (c) {
                c->addElt(n);
                n->_comp = c;
            }

            // set index
            n->set_i((size_t) i);

            // add to network
            net->addNode(n);
        }

        // remove empty compartments
        net->elideEmptyComps();

        // resize compartments to enclose contents
        //NOTE: nodes do not yet have positions, can't do this
        //net->resizeCompartments();
        net->autosizeComps();

        // add connections
        for (int i_rxn = 0; i_rxn < mod.getNumReactions(); ++i_rxn) {
            const libsbml::Reaction *rxn = mod.getReaction(i_rxn);
            auto *r = new Reaction();

            r->setId(rxn->getId());

            if (!rxn)
                throw std::logic_error("Failed to get reaction");

            // associate compartment (if one exists)
            sbnw::Compartment *c = net->findCompById(rxn->getCompartment());
            if (c) {
                c->addElt(r);
            }

            // get reactants
            for (int i_spc = 0; i_spc < rxn->getNumReactants(); ++i_spc) {
                //get the reference
                const libsbml::SpeciesReference *spc = rxn->getReactant(i_spc);
                Node *src = net->findNodeById(spc->getSpecies());
                if (!src)
                    throw std::logic_error("Invalid species reference");
                r->addSpeciesRef(src, RXN_ROLE_SUBSTRATE);
            }

            // get products
            for (int i_spc = 0; i_spc < rxn->getNumProducts(); ++i_spc) {
                //get the reference
                const libsbml::SpeciesReference *spc = rxn->getProduct(i_spc);
                Node *src = net->findNodeById(spc->getSpecies());
                if (!src)
                    throw std::logic_error("Invalid species reference");
                r->addSpeciesRef(src, RXN_ROLE_PRODUCT);
            }

            // get modifiers
            for (int i_spc = 0; i_spc < rxn->getNumModifiers(); ++i_spc) {
                //get the reference
                const libsbml::ModifierSpeciesReference *spc = rxn->getModifier(i_spc);
                Node *src = net->findNodeById(spc->getSpecies());
                if (!src)
                    throw std::logic_error("Invalid species reference");
                r->addSpeciesRef(src, RXN_ROLE_MODIFIER);
            }

            net->addReaction(r);
        }

        return net;
    }

}
