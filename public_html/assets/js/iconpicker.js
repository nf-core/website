/*!
 * Font Awesome Icon Picker
 * https://farbelous.github.io/fontawesome-iconpicker/
 *
 * @author Javi Aguilar, itsjavi.com
 * @license MIT License
 * @see https://github.com/farbelous/fontawesome-iconpicker/blob/master/LICENSE
 */


(function(e) {
    if (typeof define === "function" && define.amd) {
        define([ "jquery" ], e);
    } else {
        e(jQuery);
    }
})(function(j) {
    j.ui = j.ui || {};
    var e = j.ui.version = "1.12.1";
    (function() {
        var r, y = Math.max, x = Math.abs, s = /left|center|right/, i = /top|center|bottom/, f = /[\+\-]\d+(\.[\d]+)?%?/, l = /^\w+/, c = /%$/, a = j.fn.pos;
        function q(e, a, t) {
            return [ parseFloat(e[0]) * (c.test(e[0]) ? a / 100 : 1), parseFloat(e[1]) * (c.test(e[1]) ? t / 100 : 1) ];
        }
        function C(e, a) {
            return parseInt(j.css(e, a), 10) || 0;
        }
        function t(e) {
            var a = e[0];
            if (a.nodeType === 9) {
                return {
                    width: e.width(),
                    height: e.height(),
                    offset: {
                        top: 0,
                        left: 0
                    }
                };
            }
            if (j.isWindow(a)) {
                return {
                    width: e.width(),
                    height: e.height(),
                    offset: {
                        top: e.scrollTop(),
                        left: e.scrollLeft()
                    }
                };
            }
            if (a.preventDefault) {
                return {
                    width: 0,
                    height: 0,
                    offset: {
                        top: a.pageY,
                        left: a.pageX
                    }
                };
            }
            return {
                width: e.outerWidth(),
                height: e.outerHeight(),
                offset: e.offset()
            };
        }
        j.pos = {
            scrollbarWidth: function() {
                if (r !== undefined) {
                    return r;
                }
                var e, a, t = j("<div " + "style='display:block;position:absolute;width:50px;height:50px;overflow:hidden;'>" + "<div style='height:100px;width:auto;'></div></div>"), s = t.children()[0];
                j("body").append(t);
                e = s.offsetWidth;
                t.css("overflow", "scroll");
                a = s.offsetWidth;
                if (e === a) {
                    a = t[0].clientWidth;
                }
                t.remove();
                return r = e - a;
            },
            getScrollInfo: function(e) {
                var a = e.isWindow || e.isDocument ? "" : e.element.css("overflow-x"), t = e.isWindow || e.isDocument ? "" : e.element.css("overflow-y"), s = a === "scroll" || a === "auto" && e.width < e.element[0].scrollWidth, r = t === "scroll" || t === "auto" && e.height < e.element[0].scrollHeight;
                return {
                    width: r ? j.pos.scrollbarWidth() : 0,
                    height: s ? j.pos.scrollbarWidth() : 0
                };
            },
            getWithinInfo: function(e) {
                var a = j(e || window), t = j.isWindow(a[0]), s = !!a[0] && a[0].nodeType === 9, r = !t && !s;
                return {
                    element: a,
                    isWindow: t,
                    isDocument: s,
                    offset: r ? j(e).offset() : {
                        left: 0,
                        top: 0
                    },
                    scrollLeft: a.scrollLeft(),
                    scrollTop: a.scrollTop(),
                    width: a.outerWidth(),
                    height: a.outerHeight()
                };
            }
        };
        j.fn.pos = function(h) {
            if (!h || !h.of) {
                return a.apply(this, arguments);
            }
            h = j.extend({}, h);
            var m, p, d, u, T, e, g = j(h.of), b = j.pos.getWithinInfo(h.within), k = j.pos.getScrollInfo(b), w = (h.collision || "flip").split(" "), v = {};
            e = t(g);
            if (g[0].preventDefault) {
                h.at = "left top";
            }
            p = e.width;
            d = e.height;
            u = e.offset;
            T = j.extend({}, u);
            j.each([ "my", "at" ], function() {
                var e = (h[this] || "").split(" "), a, t;
                if (e.length === 1) {
                    e = s.test(e[0]) ? e.concat([ "center" ]) : i.test(e[0]) ? [ "center" ].concat(e) : [ "center", "center" ];
                }
                e[0] = s.test(e[0]) ? e[0] : "center";
                e[1] = i.test(e[1]) ? e[1] : "center";
                a = f.exec(e[0]);
                t = f.exec(e[1]);
                v[this] = [ a ? a[0] : 0, t ? t[0] : 0 ];
                h[this] = [ l.exec(e[0])[0], l.exec(e[1])[0] ];
            });
            if (w.length === 1) {
                w[1] = w[0];
            }
            if (h.at[0] === "right") {
                T.left += p;
            } else if (h.at[0] === "center") {
                T.left += p / 2;
            }
            if (h.at[1] === "bottom") {
                T.top += d;
            } else if (h.at[1] === "center") {
                T.top += d / 2;
            }
            m = q(v.at, p, d);
            T.left += m[0];
            T.top += m[1];
            return this.each(function() {
                var t, e, f = j(this), l = f.outerWidth(), c = f.outerHeight(), a = C(this, "marginLeft"), s = C(this, "marginTop"), r = l + a + C(this, "marginRight") + k.width, i = c + s + C(this, "marginBottom") + k.height, o = j.extend({}, T), n = q(v.my, f.outerWidth(), f.outerHeight());
                if (h.my[0] === "right") {
                    o.left -= l;
                } else if (h.my[0] === "center") {
                    o.left -= l / 2;
                }
                if (h.my[1] === "bottom") {
                    o.top -= c;
                } else if (h.my[1] === "center") {
                    o.top -= c / 2;
                }
                o.left += n[0];
                o.top += n[1];
                t = {
                    marginLeft: a,
                    marginTop: s
                };
                j.each([ "left", "top" ], function(e, a) {
                    if (j.ui.pos[w[e]]) {
                        j.ui.pos[w[e]][a](o, {
                            targetWidth: p,
                            targetHeight: d,
                            elemWidth: l,
                            elemHeight: c,
                            collisionPosition: t,
                            collisionWidth: r,
                            collisionHeight: i,
                            offset: [ m[0] + n[0], m[1] + n[1] ],
                            my: h.my,
                            at: h.at,
                            within: b,
                            elem: f
                        });
                    }
                });
                if (h.using) {
                    e = function(e) {
                        var a = u.left - o.left, t = a + p - l, s = u.top - o.top, r = s + d - c, i = {
                            target: {
                                element: g,
                                left: u.left,
                                top: u.top,
                                width: p,
                                height: d
                            },
                            element: {
                                element: f,
                                left: o.left,
                                top: o.top,
                                width: l,
                                height: c
                            },
                            horizontal: t < 0 ? "left" : a > 0 ? "right" : "center",
                            vertical: r < 0 ? "top" : s > 0 ? "bottom" : "middle"
                        };
                        if (p < l && x(a + t) < p) {
                            i.horizontal = "center";
                        }
                        if (d < c && x(s + r) < d) {
                            i.vertical = "middle";
                        }
                        if (y(x(a), x(t)) > y(x(s), x(r))) {
                            i.important = "horizontal";
                        } else {
                            i.important = "vertical";
                        }
                        h.using.call(this, e, i);
                    };
                }
                f.offset(j.extend(o, {
                    using: e
                }));
            });
        };
        j.ui.pos = {
            _trigger: function(e, a, t, s) {
                if (a.elem) {
                    a.elem.trigger({
                        type: t,
                        position: e,
                        positionData: a,
                        triggered: s
                    });
                }
            },
            fit: {
                left: function(e, a) {
                    j.ui.pos._trigger(e, a, "posCollide", "fitLeft");
                    var t = a.within, s = t.isWindow ? t.scrollLeft : t.offset.left, r = t.width, i = e.left - a.collisionPosition.marginLeft, f = s - i, l = i + a.collisionWidth - r - s, c;
                    if (a.collisionWidth > r) {
                        if (f > 0 && l <= 0) {
                            c = e.left + f + a.collisionWidth - r - s;
                            e.left += f - c;
                        } else if (l > 0 && f <= 0) {
                            e.left = s;
                        } else {
                            if (f > l) {
                                e.left = s + r - a.collisionWidth;
                            } else {
                                e.left = s;
                            }
                        }
                    } else if (f > 0) {
                        e.left += f;
                    } else if (l > 0) {
                        e.left -= l;
                    } else {
                        e.left = y(e.left - i, e.left);
                    }
                    j.ui.pos._trigger(e, a, "posCollided", "fitLeft");
                },
                top: function(e, a) {
                    j.ui.pos._trigger(e, a, "posCollide", "fitTop");
                    var t = a.within, s = t.isWindow ? t.scrollTop : t.offset.top, r = a.within.height, i = e.top - a.collisionPosition.marginTop, f = s - i, l = i + a.collisionHeight - r - s, c;
                    if (a.collisionHeight > r) {
                        if (f > 0 && l <= 0) {
                            c = e.top + f + a.collisionHeight - r - s;
                            e.top += f - c;
                        } else if (l > 0 && f <= 0) {
                            e.top = s;
                        } else {
                            if (f > l) {
                                e.top = s + r - a.collisionHeight;
                            } else {
                                e.top = s;
                            }
                        }
                    } else if (f > 0) {
                        e.top += f;
                    } else if (l > 0) {
                        e.top -= l;
                    } else {
                        e.top = y(e.top - i, e.top);
                    }
                    j.ui.pos._trigger(e, a, "posCollided", "fitTop");
                }
            },
            flip: {
                left: function(e, a) {
                    j.ui.pos._trigger(e, a, "posCollide", "flipLeft");
                    var t = a.within, s = t.offset.left + t.scrollLeft, r = t.width, i = t.isWindow ? t.scrollLeft : t.offset.left, f = e.left - a.collisionPosition.marginLeft, l = f - i, c = f + a.collisionWidth - r - i, o = a.my[0] === "left" ? -a.elemWidth : a.my[0] === "right" ? a.elemWidth : 0, n = a.at[0] === "left" ? a.targetWidth : a.at[0] === "right" ? -a.targetWidth : 0, h = -2 * a.offset[0], m, p;
                    if (l < 0) {
                        m = e.left + o + n + h + a.collisionWidth - r - s;
                        if (m < 0 || m < x(l)) {
                            e.left += o + n + h;
                        }
                    } else if (c > 0) {
                        p = e.left - a.collisionPosition.marginLeft + o + n + h - i;
                        if (p > 0 || x(p) < c) {
                            e.left += o + n + h;
                        }
                    }
                    j.ui.pos._trigger(e, a, "posCollided", "flipLeft");
                },
                top: function(e, a) {
                    j.ui.pos._trigger(e, a, "posCollide", "flipTop");
                    var t = a.within, s = t.offset.top + t.scrollTop, r = t.height, i = t.isWindow ? t.scrollTop : t.offset.top, f = e.top - a.collisionPosition.marginTop, l = f - i, c = f + a.collisionHeight - r - i, o = a.my[1] === "top", n = o ? -a.elemHeight : a.my[1] === "bottom" ? a.elemHeight : 0, h = a.at[1] === "top" ? a.targetHeight : a.at[1] === "bottom" ? -a.targetHeight : 0, m = -2 * a.offset[1], p, d;
                    if (l < 0) {
                        d = e.top + n + h + m + a.collisionHeight - r - s;
                        if (d < 0 || d < x(l)) {
                            e.top += n + h + m;
                        }
                    } else if (c > 0) {
                        p = e.top - a.collisionPosition.marginTop + n + h + m - i;
                        if (p > 0 || x(p) < c) {
                            e.top += n + h + m;
                        }
                    }
                    j.ui.pos._trigger(e, a, "posCollided", "flipTop");
                }
            },
            flipfit: {
                left: function() {
                    j.ui.pos.flip.left.apply(this, arguments);
                    j.ui.pos.fit.left.apply(this, arguments);
                },
                top: function() {
                    j.ui.pos.flip.top.apply(this, arguments);
                    j.ui.pos.fit.top.apply(this, arguments);
                }
            }
        };
        (function() {
            var e, a, t, s, r, i = document.getElementsByTagName("body")[0], f = document.createElement("div");
            e = document.createElement(i ? "div" : "body");
            t = {
                visibility: "hidden",
                width: 0,
                height: 0,
                border: 0,
                margin: 0,
                background: "none"
            };
            if (i) {
                j.extend(t, {
                    position: "absolute",
                    left: "-1000px",
                    top: "-1000px"
                });
            }
            for (r in t) {
                e.style[r] = t[r];
            }
            e.appendChild(f);
            a = i || document.documentElement;
            a.insertBefore(e, a.firstChild);
            f.style.cssText = "position: absolute; left: 10.7432222px;";
            s = j(f).offset().left;
            j.support.offsetFractions = s > 10 && s < 11;
            e.innerHTML = "";
            a.removeChild(e);
        })();
    })();
    var a = j.ui.position;
});

(function(e) {
    "use strict";
    if (typeof define === "function" && define.amd) {
        define([ "jquery" ], e);
    } else if (window.jQuery && !window.jQuery.fn.iconpicker) {
        e(window.jQuery);
    }
})(function(c) {
    "use strict";
    var f = {
        isEmpty: function(e) {
            return e === false || e === "" || e === null || e === undefined;
        },
        isEmptyObject: function(e) {
            return this.isEmpty(e) === true || e.length === 0;
        },
        isElement: function(e) {
            return c(e).length > 0;
        },
        isString: function(e) {
            return typeof e === "string" || e instanceof String;
        },
        isArray: function(e) {
            return c.isArray(e);
        },
        inArray: function(e, a) {
            return c.inArray(e, a) !== -1;
        },
        throwError: function(e) {
            throw "Font Awesome Icon Picker Exception: " + e;
        }
    };
    var t = function(e, a) {
        this._id = t._idCounter++;
        this.element = c(e).addClass("iconpicker-element");
        this._trigger("iconpickerCreate", {
            iconpickerValue: this.iconpickerValue
        });
        this.options = c.extend({}, t.defaultOptions, this.element.data(), a);
        this.options.templates = c.extend({}, t.defaultOptions.templates, this.options.templates);
        this.options.originalPlacement = this.options.placement;
        this.container = f.isElement(this.options.container) ? c(this.options.container) : false;
        if (this.container === false) {
            if (this.element.is(".dropdown-toggle")) {
                this.container = c("~ .dropdown-menu:first", this.element);
            } else {
                this.container = this.element.is("input,textarea,button,.btn") ? this.element.parent() : this.element;
            }
        }
        this.container.addClass("iconpicker-container");
        if (this.isDropdownMenu()) {
            this.options.placement = "inline";
        }
        this.input = this.element.is("input,textarea") ? this.element.addClass("iconpicker-input") : false;
        if (this.input === false) {
            this.input = this.container.find(this.options.input);
            if (!this.input.is("input,textarea")) {
                this.input = false;
            }
        }
        this.component = this.isDropdownMenu() ? this.container.parent().find(this.options.component) : this.container.find(this.options.component);
        if (this.component.length === 0) {
            this.component = false;
        } else {
            this.component.find("i").addClass("iconpicker-component");
        }
        this._createPopover();
        this._createIconpicker();
        if (this.getAcceptButton().length === 0) {
            this.options.mustAccept = false;
        }
        if (this.isInputGroup()) {
            this.container.parent().append(this.popover);
        } else {
            this.container.append(this.popover);
        }
        this._bindElementEvents();
        this._bindWindowEvents();
        this.update(this.options.selected);
        if (this.isInline()) {
            this.show();
        }
        this._trigger("iconpickerCreated", {
            iconpickerValue: this.iconpickerValue
        });
    };
    t._idCounter = 0;
    t.defaultOptions = {
        title: false,
        selected: false,
        defaultValue: false,
        placement: "bottom",
        collision: "none",
        animation: true,
        hideOnSelect: false,
        showFooter: false,
        searchInFooter: false,
        mustAccept: false,
        selectedCustomClass: "bg-primary",
        icons: [],
        fullClassFormatter: function(e) {
            return e;
        },
        input: "input,.iconpicker-input",
        inputSearch: false,
        container: false,
        component: ".input-group-addon,.iconpicker-component",
        templates: {
            popover: '<div class="iconpicker-popover popover"><div class="arrow"></div>' + '<div class="popover-title"></div><div class="popover-content"></div></div>',
            footer: '<div class="popover-footer"></div>',
            buttons: '<button class="iconpicker-btn iconpicker-btn-cancel btn btn-default btn-sm">Cancel</button>' + ' <button class="iconpicker-btn iconpicker-btn-accept btn btn-primary btn-sm">Accept</button>',
            search: '<input type="search" class="form-control iconpicker-search" placeholder="Type to filter" />',
            iconpicker: '<div class="iconpicker"><div class="iconpicker-items"></div></div>',
            iconpickerItem: '<a role="button" href="javascript:;" class="iconpicker-item"><i></i></a>'
        }
    };
    t.batch = function(e, a) {
        var t = Array.prototype.slice.call(arguments, 2);
        return c(e).each(function() {
            var e = c(this).data("iconpicker");
            if (!!e) {
                e[a].apply(e, t);
            }
        });
    };
    t.prototype = {
        constructor: t,
        options: {},
        _id: 0,
        _trigger: function(e, a) {
            a = a || {};
            this.element.trigger(c.extend({
                type: e,
                iconpickerInstance: this
            }, a));
        },
        _createPopover: function() {
            this.popover = c(this.options.templates.popover);
            var e = this.popover.find(".popover-title");
            if (!!this.options.title) {
                e.append(c('<div class="popover-title-text">' + this.options.title + "</div>"));
            }
            if (this.hasSeparatedSearchInput() && !this.options.searchInFooter) {
                e.append(this.options.templates.search);
            } else if (!this.options.title) {
                e.remove();
            }
            if (this.options.showFooter && !f.isEmpty(this.options.templates.footer)) {
                var a = c(this.options.templates.footer);
                if (this.hasSeparatedSearchInput() && this.options.searchInFooter) {
                    a.append(c(this.options.templates.search));
                }
                if (!f.isEmpty(this.options.templates.buttons)) {
                    a.append(c(this.options.templates.buttons));
                }
                this.popover.append(a);
            }
            if (this.options.animation === true) {
                this.popover.addClass("fade");
            }
            return this.popover;
        },
        _createIconpicker: function() {
            var t = this;
            this.iconpicker = c(this.options.templates.iconpicker);
            var e = function(e) {
                var a = c(this);
                if (a.is("i")) {
                    a = a.parent();
                }
                t._trigger("iconpickerSelect", {
                    iconpickerItem: a,
                    iconpickerValue: t.iconpickerValue
                });
                if (t.options.mustAccept === false) {
                    t.update(a.data("iconpickerValue"));
                    t._trigger("iconpickerSelected", {
                        iconpickerItem: this,
                        iconpickerValue: t.iconpickerValue
                    });
                } else {
                    t.update(a.data("iconpickerValue"), true);
                }
                if (t.options.hideOnSelect && t.options.mustAccept === false) {
                    t.hide();
                }
            };
            var a = c(this.options.templates.iconpickerItem);
            var s = [];
            for (var r in this.options.icons) {
                if (typeof this.options.icons[r].title === "string") {
                    var i = a.clone();
                    i.find("i").addClass(this.options.fullClassFormatter(this.options.icons[r].title));
                    i.data("iconpickerValue", this.options.icons[r].title).on("click.iconpicker", e);
                    i.attr("title", "." + this.options.icons[r].title);
                    if (this.options.icons[r].searchTerms.length > 0) {
                        var f = "";
                        for (var l = 0; l < this.options.icons[r].searchTerms.length; l++) {
                            f = f + this.options.icons[r].searchTerms[l] + " ";
                        }
                        i.attr("data-search-terms", f);
                    }
                    s.push(i);
                }
            }
            this.iconpicker.find(".iconpicker-items").append(s);
            this.popover.find(".popover-content").append(this.iconpicker);
            return this.iconpicker;
        },
        _isEventInsideIconpicker: function(e) {
            var a = c(e.target);
            if ((!a.hasClass("iconpicker-element") || a.hasClass("iconpicker-element") && !a.is(this.element)) && a.parents(".iconpicker-popover").length === 0) {
                return false;
            }
            return true;
        },
        _bindElementEvents: function() {
            var a = this;
            this.getSearchInput().on("keyup.iconpicker", function() {
                a.filter(c(this).val().toLowerCase());
            });
            this.getAcceptButton().on("click.iconpicker", function() {
                var e = a.iconpicker.find(".iconpicker-selected").get(0);
                a.update(a.iconpickerValue);
                a._trigger("iconpickerSelected", {
                    iconpickerItem: e,
                    iconpickerValue: a.iconpickerValue
                });
                if (!a.isInline()) {
                    a.hide();
                }
            });
            this.getCancelButton().on("click.iconpicker", function() {
                if (!a.isInline()) {
                    a.hide();
                }
            });
            this.element.on("focus.iconpicker", function(e) {
                a.show();
                e.stopPropagation();
            });
            if (this.hasComponent()) {
                this.component.on("click.iconpicker", function() {
                    a.toggle();
                });
            }
            if (this.hasInput()) {
                this.input.on("keyup.iconpicker", function(e) {
                    if (!f.inArray(e.keyCode, [ 38, 40, 37, 39, 16, 17, 18, 9, 8, 91, 93, 20, 46, 186, 190, 46, 78, 188, 44, 86 ])) {
                        a.update();
                    } else {
                        a._updateFormGroupStatus(a.getValid(this.value) !== false);
                    }
                    if (a.options.inputSearch === true) {
                        a.filter(c(this).val().toLowerCase());
                    }
                });
            }
        },
        _bindWindowEvents: function() {
            var e = c(window.document);
            var a = this;
            var t = ".iconpicker.inst" + this._id;
            c(window).on("resize.iconpicker" + t + " orientationchange.iconpicker" + t, function(e) {
                if (a.popover.hasClass("in")) {
                    a.updatePlacement();
                }
            });
            if (!a.isInline()) {
                e.on("mouseup" + t, function(e) {
                    if (!a._isEventInsideIconpicker(e) && !a.isInline()) {
                        a.hide();
                    }
                });
            }
        },
        _unbindElementEvents: function() {
            this.popover.off(".iconpicker");
            this.element.off(".iconpicker");
            if (this.hasInput()) {
                this.input.off(".iconpicker");
            }
            if (this.hasComponent()) {
                this.component.off(".iconpicker");
            }
            if (this.hasContainer()) {
                this.container.off(".iconpicker");
            }
        },
        _unbindWindowEvents: function() {
            c(window).off(".iconpicker.inst" + this._id);
            c(window.document).off(".iconpicker.inst" + this._id);
        },
        updatePlacement: function(e, a) {
            e = e || this.options.placement;
            this.options.placement = e;
            a = a || this.options.collision;
            a = a === true ? "flip" : a;
            var t = {
                at: "right bottom",
                my: "right top",
                of: this.hasInput() && !this.isInputGroup() ? this.input : this.container,
                collision: a === true ? "flip" : a,
                within: window
            };
            this.popover.removeClass("inline topLeftCorner topLeft top topRight topRightCorner " + "rightTop right rightBottom bottomRight bottomRightCorner " + "bottom bottomLeft bottomLeftCorner leftBottom left leftTop");
            if (typeof e === "object") {
                return this.popover.pos(c.extend({}, t, e));
            }
            switch (e) {
              case "inline":
                {
                    t = false;
                }
                break;

              case "topLeftCorner":
                {
                    t.my = "right bottom";
                    t.at = "left top";
                }
                break;

              case "topLeft":
                {
                    t.my = "left bottom";
                    t.at = "left top";
                }
                break;

              case "top":
                {
                    t.my = "center bottom";
                    t.at = "center top";
                }
                break;

              case "topRight":
                {
                    t.my = "right bottom";
                    t.at = "right top";
                }
                break;

              case "topRightCorner":
                {
                    t.my = "left bottom";
                    t.at = "right top";
                }
                break;

              case "rightTop":
                {
                    t.my = "left bottom";
                    t.at = "right center";
                }
                break;

              case "right":
                {
                    t.my = "left center";
                    t.at = "right center";
                }
                break;

              case "rightBottom":
                {
                    t.my = "left top";
                    t.at = "right center";
                }
                break;

              case "bottomRightCorner":
                {
                    t.my = "left top";
                    t.at = "right bottom";
                }
                break;

              case "bottomRight":
                {
                    t.my = "right top";
                    t.at = "right bottom";
                }
                break;

              case "bottom":
                {
                    t.my = "center top";
                    t.at = "center bottom";
                }
                break;

              case "bottomLeft":
                {
                    t.my = "left top";
                    t.at = "left bottom";
                }
                break;

              case "bottomLeftCorner":
                {
                    t.my = "right top";
                    t.at = "left bottom";
                }
                break;

              case "leftBottom":
                {
                    t.my = "right top";
                    t.at = "left center";
                }
                break;

              case "left":
                {
                    t.my = "right center";
                    t.at = "left center";
                }
                break;

              case "leftTop":
                {
                    t.my = "right bottom";
                    t.at = "left center";
                }
                break;

              default:
                {
                    return false;
                }
                break;
            }
            this.popover.css({
                display: this.options.placement === "inline" ? "" : "block"
            });
            if (t !== false) {
                this.popover.pos(t).css("maxWidth", c(window).width() - this.container.offset().left - 5);
            } else {
                this.popover.css({
                    top: "auto",
                    right: "auto",
                    bottom: "auto",
                    left: "auto",
                    maxWidth: "none"
                });
            }
            this.popover.addClass(this.options.placement);
            return true;
        },
        _updateComponents: function() {
            this.iconpicker.find(".iconpicker-item.iconpicker-selected").removeClass("iconpicker-selected " + this.options.selectedCustomClass);
            if (this.iconpickerValue) {
                this.iconpicker.find("." + this.options.fullClassFormatter(this.iconpickerValue).replace(/ /g, ".")).parent().addClass("iconpicker-selected " + this.options.selectedCustomClass);
            }
            if (this.hasComponent()) {
                var e = this.component.find("i");
                if (e.length > 0) {
                    e.attr("class", this.options.fullClassFormatter(this.iconpickerValue));
                } else {
                    this.component.html(this.getHtml());
                }
            }
        },
        _updateFormGroupStatus: function(e) {
            if (this.hasInput()) {
                if (e !== false) {
                    this.input.parents(".form-group:first").removeClass("has-error");
                } else {
                    this.input.parents(".form-group:first").addClass("has-error");
                }
                return true;
            }
            return false;
        },
        getValid: function(e) {
            if (!f.isString(e)) {
                e = "";
            }
            var a = e === "";
            e = c.trim(e);
            var t = false;
            for (var s = 0; s < this.options.icons.length; s++) {
                if (this.options.icons[s].title === e) {
                    t = true;
                    break;
                }
            }
            if (t || a) {
                return e;
            }
            return false;
        },
        setValue: function(e) {
            var a = this.getValid(e);
            if (a !== false) {
                this.iconpickerValue = a;
                this._trigger("iconpickerSetValue", {
                    iconpickerValue: a
                });
                return this.iconpickerValue;
            } else {
                this._trigger("iconpickerInvalid", {
                    iconpickerValue: e
                });
                return false;
            }
        },
        getHtml: function() {
            return '<i class="' + this.options.fullClassFormatter(this.iconpickerValue) + '"></i>';
        },
        setSourceValue: function(e) {
            e = this.setValue(e);
            if (e !== false && e !== "") {
                if (this.hasInput()) {
                    this.input.val(this.iconpickerValue);
                } else {
                    this.element.data("iconpickerValue", this.iconpickerValue);
                }
                this._trigger("iconpickerSetSourceValue", {
                    iconpickerValue: e
                });
            }
            return e;
        },
        getSourceValue: function(e) {
            e = e || this.options.defaultValue;
            var a = e;
            if (this.hasInput()) {
                a = this.input.val();
            } else {
                a = this.element.data("iconpickerValue");
            }
            if (a === undefined || a === "" || a === null || a === false) {
                a = e;
            }
            return a;
        },
        hasInput: function() {
            return this.input !== false;
        },
        isInputSearch: function() {
            return this.hasInput() && this.options.inputSearch === true;
        },
        isInputGroup: function() {
            return this.container.is(".input-group");
        },
        isDropdownMenu: function() {
            return this.container.is(".dropdown-menu");
        },
        hasSeparatedSearchInput: function() {
            return this.options.templates.search !== false && !this.isInputSearch();
        },
        hasComponent: function() {
            return this.component !== false;
        },
        hasContainer: function() {
            return this.container !== false;
        },
        getAcceptButton: function() {
            return this.popover.find(".iconpicker-btn-accept");
        },
        getCancelButton: function() {
            return this.popover.find(".iconpicker-btn-cancel");
        },
        getSearchInput: function() {
            return this.popover.find(".iconpicker-search");
        },
        filter: function(r) {
            if (f.isEmpty(r)) {
                this.iconpicker.find(".iconpicker-item").show();
                return c(false);
            } else {
                var i = [];
                this.iconpicker.find(".iconpicker-item").each(function() {
                    var e = c(this);
                    var a = e.attr("title").toLowerCase();
                    var t = e.attr("data-search-terms") ? e.attr("data-search-terms").toLowerCase() : "";
                    a = a + " " + t;
                    var s = false;
                    try {
                        s = new RegExp("(^|\\W)" + r, "g");
                    } catch (e) {
                        s = false;
                    }
                    if (s !== false && a.match(s)) {
                        i.push(e);
                        e.show();
                    } else {
                        e.hide();
                    }
                });
                return i;
            }
        },
        show: function() {
            if (this.popover.hasClass("in")) {
                return false;
            }
            c.iconpicker.batch(c(".iconpicker-popover.in:not(.inline)").not(this.popover), "hide");
            this._trigger("iconpickerShow", {
                iconpickerValue: this.iconpickerValue
            });
            this.updatePlacement();
            this.popover.addClass("in");
            setTimeout(c.proxy(function() {
                this.popover.css("display", this.isInline() ? "" : "block");
                this._trigger("iconpickerShown", {
                    iconpickerValue: this.iconpickerValue
                });
            }, this), this.options.animation ? 300 : 1);
        },
        hide: function() {
            if (!this.popover.hasClass("in")) {
                return false;
            }
            this._trigger("iconpickerHide", {
                iconpickerValue: this.iconpickerValue
            });
            this.popover.removeClass("in");
            setTimeout(c.proxy(function() {
                this.popover.css("display", "none");
                this.getSearchInput().val("");
                this.filter("");
                this._trigger("iconpickerHidden", {
                    iconpickerValue: this.iconpickerValue
                });
            }, this), this.options.animation ? 300 : 1);
        },
        toggle: function() {
            if (this.popover.is(":visible")) {
                this.hide();
            } else {
                this.show(true);
            }
        },
        update: function(e, a) {
            e = e ? e : this.getSourceValue(this.iconpickerValue);
            this._trigger("iconpickerUpdate", {
                iconpickerValue: this.iconpickerValue
            });
            if (a === true) {
                e = this.setValue(e);
            } else {
                e = this.setSourceValue(e);
                this._updateFormGroupStatus(e !== false);
            }
            if (e !== false) {
                this._updateComponents();
            }
            this._trigger("iconpickerUpdated", {
                iconpickerValue: this.iconpickerValue
            });
            return e;
        },
        destroy: function() {
            this._trigger("iconpickerDestroy", {
                iconpickerValue: this.iconpickerValue
            });
            this.element.removeData("iconpicker").removeData("iconpickerValue").removeClass("iconpicker-element");
            this._unbindElementEvents();
            this._unbindWindowEvents();
            c(this.popover).remove();
            this._trigger("iconpickerDestroyed", {
                iconpickerValue: this.iconpickerValue
            });
        },
        disable: function() {
            if (this.hasInput()) {
                this.input.prop("disabled", true);
                return true;
            }
            return false;
        },
        enable: function() {
            if (this.hasInput()) {
                this.input.prop("disabled", false);
                return true;
            }
            return false;
        },
        isDisabled: function() {
            if (this.hasInput()) {
                return this.input.prop("disabled") === true;
            }
            return false;
        },
        isInline: function() {
            return this.options.placement === "inline" || this.popover.hasClass("inline");
        }
    };
    c.iconpicker = t;
    c.fn.iconpicker = function(a) {
        return this.each(function() {
            var e = c(this);
            if (!e.data("iconpicker")) {
                e.data("iconpicker", new t(this, typeof a === "object" ? a : {}));
            }
        });
    };
    t.defaultOptions = c.extend(t.defaultOptions, {"icons": [{"title": "fab fa-500px", "searchTerms": []}, {"title": "fab fa-accessible-icon", "searchTerms": ["accessibility", "handicap", "person", "wheelchair", "wheelchair-alt"]}, {"title": "fab fa-accusoft", "searchTerms": []}, {"title": "fab fa-acquisitions-incorporated", "searchTerms": ["Dungeons & Dragons", "d&d", "dnd", "fantasy", "game", "gaming", "tabletop"]}, {"title": "fas fa-ad", "searchTerms": ["advertisement", "media", "newspaper", "promotion", "publicity"]}, {"title": "fas fa-address-book", "searchTerms": ["contact", "directory", "index", "little black book", "rolodex"]}, {"title": "far fa-address-book", "searchTerms": ["contact", "directory", "index", "little black book", "rolodex"]}, {"title": "fas fa-address-card", "searchTerms": ["about", "contact", "id", "identification", "postcard", "profile"]}, {"title": "far fa-address-card", "searchTerms": ["about", "contact", "id", "identification", "postcard", "profile"]}, {"title": "fas fa-adjust", "searchTerms": ["contrast", "dark", "light", "saturation"]}, {"title": "fab fa-adn", "searchTerms": []}, {"title": "fab fa-adobe", "searchTerms": ["acrobat", "app", "design", "illustrator", "indesign", "photoshop"]}, {"title": "fab fa-adversal", "searchTerms": []}, {"title": "fab fa-affiliatetheme", "searchTerms": []}, {"title": "fas fa-air-freshener", "searchTerms": ["car", "deodorize", "fresh", "pine", "scent"]}, {"title": "fab fa-airbnb", "searchTerms": []}, {"title": "fab fa-algolia", "searchTerms": []}, {"title": "fas fa-align-center", "searchTerms": ["format", "middle", "paragraph", "text"]}, {"title": "fas fa-align-justify", "searchTerms": ["format", "paragraph", "text"]}, {"title": "fas fa-align-left", "searchTerms": ["format", "paragraph", "text"]}, {"title": "fas fa-align-right", "searchTerms": ["format", "paragraph", "text"]}, {"title": "fab fa-alipay", "searchTerms": []}, {"title": "fas fa-allergies", "searchTerms": ["allergy", "freckles", "hand", "hives", "pox", "skin", "spots"]}, {"title": "fab fa-amazon", "searchTerms": []}, {"title": "fab fa-amazon-pay", "searchTerms": []}, {"title": "fas fa-ambulance", "searchTerms": ["emergency", "emt", "er", "help", "hospital", "support", "vehicle"]}, {"title": "fas fa-american-sign-language-interpreting", "searchTerms": ["asl", "deaf", "finger", "hand", "interpret", "speak"]}, {"title": "fab fa-amilia", "searchTerms": []}, {"title": "fas fa-anchor", "searchTerms": ["berth", "boat", "dock", "embed", "link", "maritime", "moor", "secure"]}, {"title": "fab fa-android", "searchTerms": ["robot"]}, {"title": "fab fa-angellist", "searchTerms": []}, {"title": "fas fa-angle-double-down", "searchTerms": ["arrows", "caret", "download", "expand"]}, {"title": "fas fa-angle-double-left", "searchTerms": ["arrows", "back", "caret", "laquo", "previous", "quote"]}, {"title": "fas fa-angle-double-right", "searchTerms": ["arrows", "caret", "forward", "more", "next", "quote", "raquo"]}, {"title": "fas fa-angle-double-up", "searchTerms": ["arrows", "caret", "collapse", "upload"]}, {"title": "fas fa-angle-down", "searchTerms": ["arrow", "caret", "download", "expand"]}, {"title": "fas fa-angle-left", "searchTerms": ["arrow", "back", "caret", "less", "previous"]}, {"title": "fas fa-angle-right", "searchTerms": ["arrow", "care", "forward", "more", "next"]}, {"title": "fas fa-angle-up", "searchTerms": ["arrow", "caret", "collapse", "upload"]}, {"title": "fas fa-angry", "searchTerms": ["disapprove", "emoticon", "face", "mad", "upset"]}, {"title": "far fa-angry", "searchTerms": ["disapprove", "emoticon", "face", "mad", "upset"]}, {"title": "fab fa-angrycreative", "searchTerms": []}, {"title": "fab fa-angular", "searchTerms": []}, {"title": "fas fa-ankh", "searchTerms": ["amulet", "copper", "coptic christianity", "copts", "crux ansata", "egypt", "venus"]}, {"title": "fab fa-app-store", "searchTerms": []}, {"title": "fab fa-app-store-ios", "searchTerms": []}, {"title": "fab fa-apper", "searchTerms": []}, {"title": "fab fa-apple", "searchTerms": ["fruit", "ios", "mac", "operating system", "os", "osx"]}, {"title": "fas fa-apple-alt", "searchTerms": ["fall", "fruit", "fuji", "macintosh", "orchard", "seasonal", "vegan"]}, {"title": "fab fa-apple-pay", "searchTerms": []}, {"title": "fas fa-archive", "searchTerms": ["box", "package", "save", "storage"]}, {"title": "fas fa-archway", "searchTerms": ["arc", "monument", "road", "street", "tunnel"]}, {"title": "fas fa-arrow-alt-circle-down", "searchTerms": ["arrow-circle-o-down", "download"]}, {"title": "far fa-arrow-alt-circle-down", "searchTerms": ["arrow-circle-o-down", "download"]}, {"title": "fas fa-arrow-alt-circle-left", "searchTerms": ["arrow-circle-o-left", "back", "previous"]}, {"title": "far fa-arrow-alt-circle-left", "searchTerms": ["arrow-circle-o-left", "back", "previous"]}, {"title": "fas fa-arrow-alt-circle-right", "searchTerms": ["arrow-circle-o-right", "forward", "next"]}, {"title": "far fa-arrow-alt-circle-right", "searchTerms": ["arrow-circle-o-right", "forward", "next"]}, {"title": "fas fa-arrow-alt-circle-up", "searchTerms": ["arrow-circle-o-up"]}, {"title": "far fa-arrow-alt-circle-up", "searchTerms": ["arrow-circle-o-up"]}, {"title": "fas fa-arrow-circle-down", "searchTerms": ["download"]}, {"title": "fas fa-arrow-circle-left", "searchTerms": ["back", "previous"]}, {"title": "fas fa-arrow-circle-right", "searchTerms": ["forward", "next"]}, {"title": "fas fa-arrow-circle-up", "searchTerms": ["upload"]}, {"title": "fas fa-arrow-down", "searchTerms": ["download"]}, {"title": "fas fa-arrow-left", "searchTerms": ["back", "previous"]}, {"title": "fas fa-arrow-right", "searchTerms": ["forward", "next"]}, {"title": "fas fa-arrow-up", "searchTerms": ["forward", "upload"]}, {"title": "fas fa-arrows-alt", "searchTerms": ["arrow", "arrows", "bigger", "enlarge", "expand", "fullscreen", "move", "position", "reorder", "resize"]}, {"title": "fas fa-arrows-alt-h", "searchTerms": ["arrows-h", "expand", "horizontal", "landscape", "resize", "wide"]}, {"title": "fas fa-arrows-alt-v", "searchTerms": ["arrows-v", "expand", "portrait", "resize", "tall", "vertical"]}, {"title": "fab fa-artstation", "searchTerms": []}, {"title": "fas fa-assistive-listening-systems", "searchTerms": ["amplify", "audio", "deaf", "ear", "headset", "hearing", "sound"]}, {"title": "fas fa-asterisk", "searchTerms": ["annotation", "details", "reference", "star"]}, {"title": "fab fa-asymmetrik", "searchTerms": []}, {"title": "fas fa-at", "searchTerms": ["address", "author", "e-mail", "email", "handle"]}, {"title": "fas fa-atlas", "searchTerms": ["book", "directions", "geography", "globe", "map", "travel", "wayfinding"]}, {"title": "fab fa-atlassian", "searchTerms": []}, {"title": "fas fa-atom", "searchTerms": ["atheism", "chemistry", "electron", "ion", "isotope", "neutron", "nuclear", "proton", "science"]}, {"title": "fab fa-audible", "searchTerms": []}, {"title": "fas fa-audio-description", "searchTerms": ["blind", "narration", "video", "visual"]}, {"title": "fab fa-autoprefixer", "searchTerms": []}, {"title": "fab fa-avianex", "searchTerms": []}, {"title": "fab fa-aviato", "searchTerms": []}, {"title": "fas fa-award", "searchTerms": ["honor", "praise", "prize", "recognition", "ribbon", "trophy"]}, {"title": "fab fa-aws", "searchTerms": []}, {"title": "fas fa-baby", "searchTerms": ["child", "diaper", "doll", "human", "infant", "kid", "offspring", "person", "sprout"]}, {"title": "fas fa-baby-carriage", "searchTerms": ["buggy", "carrier", "infant", "push", "stroller", "transportation", "walk", "wheels"]}, {"title": "fas fa-backspace", "searchTerms": ["command", "delete", "erase", "keyboard", "undo"]}, {"title": "fas fa-backward", "searchTerms": ["previous", "rewind"]}, {"title": "fas fa-bacon", "searchTerms": ["blt", "breakfast", "ham", "lard", "meat", "pancetta", "pork", "rasher"]}, {"title": "fas fa-bahai", "searchTerms": ["bahai", "bahá'í", "star"]}, {"title": "fas fa-balance-scale", "searchTerms": ["balanced", "justice", "legal", "measure", "weight"]}, {"title": "fas fa-balance-scale-left", "searchTerms": ["justice", "legal", "measure", "unbalanced", "weight"]}, {"title": "fas fa-balance-scale-right", "searchTerms": ["justice", "legal", "measure", "unbalanced", "weight"]}, {"title": "fas fa-ban", "searchTerms": ["abort", "ban", "block", "cancel", "delete", "hide", "prohibit", "remove", "stop", "trash"]}, {"title": "fas fa-band-aid", "searchTerms": ["bandage", "boo boo", "first aid", "ouch"]}, {"title": "fab fa-bandcamp", "searchTerms": []}, {"title": "fas fa-barcode", "searchTerms": ["info", "laser", "price", "scan", "upc"]}, {"title": "fas fa-bars", "searchTerms": ["checklist", "drag", "hamburger", "list", "menu", "nav", "navigation", "ol", "reorder", "settings", "todo", "ul"]}, {"title": "fas fa-baseball-ball", "searchTerms": ["foul", "hardball", "league", "leather", "mlb", "softball", "sport"]}, {"title": "fas fa-basketball-ball", "searchTerms": ["dribble", "dunk", "hoop", "nba"]}, {"title": "fas fa-bath", "searchTerms": ["clean", "shower", "tub", "wash"]}, {"title": "fas fa-battery-empty", "searchTerms": ["charge", "dead", "power", "status"]}, {"title": "fas fa-battery-full", "searchTerms": ["charge", "power", "status"]}, {"title": "fas fa-battery-half", "searchTerms": ["charge", "power", "status"]}, {"title": "fas fa-battery-quarter", "searchTerms": ["charge", "low", "power", "status"]}, {"title": "fas fa-battery-three-quarters", "searchTerms": ["charge", "power", "status"]}, {"title": "fab fa-battle-net", "searchTerms": []}, {"title": "fas fa-bed", "searchTerms": ["lodging", "mattress", "rest", "sleep", "travel"]}, {"title": "fas fa-beer", "searchTerms": ["alcohol", "ale", "bar", "beverage", "brewery", "drink", "lager", "liquor", "mug", "stein"]}, {"title": "fab fa-behance", "searchTerms": []}, {"title": "fab fa-behance-square", "searchTerms": []}, {"title": "fas fa-bell", "searchTerms": ["alarm", "alert", "chime", "notification", "reminder"]}, {"title": "far fa-bell", "searchTerms": ["alarm", "alert", "chime", "notification", "reminder"]}, {"title": "fas fa-bell-slash", "searchTerms": ["alert", "cancel", "disabled", "notification", "off", "reminder"]}, {"title": "far fa-bell-slash", "searchTerms": ["alert", "cancel", "disabled", "notification", "off", "reminder"]}, {"title": "fas fa-bezier-curve", "searchTerms": ["curves", "illustrator", "lines", "path", "vector"]}, {"title": "fas fa-bible", "searchTerms": ["book", "catholicism", "christianity", "god", "holy"]}, {"title": "fas fa-bicycle", "searchTerms": ["bike", "gears", "pedal", "transportation", "vehicle"]}, {"title": "fas fa-biking", "searchTerms": ["bicycle", "bike", "cycle", "cycling", "ride", "wheel"]}, {"title": "fab fa-bimobject", "searchTerms": []}, {"title": "fas fa-binoculars", "searchTerms": ["glasses", "magnify", "scenic", "spyglass", "view"]}, {"title": "fas fa-biohazard", "searchTerms": ["danger", "dangerous", "hazmat", "medical", "radioactive", "toxic", "waste", "zombie"]}, {"title": "fas fa-birthday-cake", "searchTerms": ["anniversary", "bakery", "candles", "celebration", "dessert", "frosting", "holiday", "party", "pastry"]}, {"title": "fab fa-bitbucket", "searchTerms": ["atlassian", "bitbucket-square", "git"]}, {"title": "fab fa-bitcoin", "searchTerms": []}, {"title": "fab fa-bity", "searchTerms": []}, {"title": "fab fa-black-tie", "searchTerms": []}, {"title": "fab fa-blackberry", "searchTerms": []}, {"title": "fas fa-blender", "searchTerms": ["cocktail", "milkshake", "mixer", "puree", "smoothie"]}, {"title": "fas fa-blender-phone", "searchTerms": ["appliance", "cocktail", "communication", "fantasy", "milkshake", "mixer", "puree", "silly", "smoothie"]}, {"title": "fas fa-blind", "searchTerms": ["cane", "disability", "person", "sight"]}, {"title": "fas fa-blog", "searchTerms": ["journal", "log", "online", "personal", "post", "web 2.0", "wordpress", "writing"]}, {"title": "fab fa-blogger", "searchTerms": []}, {"title": "fab fa-blogger-b", "searchTerms": []}, {"title": "fab fa-bluetooth", "searchTerms": []}, {"title": "fab fa-bluetooth-b", "searchTerms": []}, {"title": "fas fa-bold", "searchTerms": ["emphasis", "format", "text"]}, {"title": "fas fa-bolt", "searchTerms": ["electricity", "lightning", "weather", "zap"]}, {"title": "fas fa-bomb", "searchTerms": ["error", "explode", "fuse", "grenade", "warning"]}, {"title": "fas fa-bone", "searchTerms": ["calcium", "dog", "skeletal", "skeleton", "tibia"]}, {"title": "fas fa-bong", "searchTerms": ["aparatus", "cannabis", "marijuana", "pipe", "smoke", "smoking"]}, {"title": "fas fa-book", "searchTerms": ["diary", "documentation", "journal", "library", "read"]}, {"title": "fas fa-book-dead", "searchTerms": ["Dungeons & Dragons", "crossbones", "d&d", "dark arts", "death", "dnd", "documentation", "evil", "fantasy", "halloween", "holiday", "necronomicon", "read", "skull", "spell"]}, {"title": "fas fa-book-medical", "searchTerms": ["diary", "documentation", "health", "history", "journal", "library", "read", "record"]}, {"title": "fas fa-book-open", "searchTerms": ["flyer", "library", "notebook", "open book", "pamphlet", "reading"]}, {"title": "fas fa-book-reader", "searchTerms": ["flyer", "library", "notebook", "open book", "pamphlet", "reading"]}, {"title": "fas fa-bookmark", "searchTerms": ["favorite", "marker", "read", "remember", "save"]}, {"title": "far fa-bookmark", "searchTerms": ["favorite", "marker", "read", "remember", "save"]}, {"title": "fab fa-bootstrap", "searchTerms": []}, {"title": "fas fa-border-all", "searchTerms": ["cell", "grid", "outline", "stroke", "table"]}, {"title": "fas fa-border-none", "searchTerms": ["cell", "grid", "outline", "stroke", "table"]}, {"title": "fas fa-border-style", "searchTerms": []}, {"title": "fas fa-bowling-ball", "searchTerms": ["alley", "candlepin", "gutter", "lane", "strike", "tenpin"]}, {"title": "fas fa-box", "searchTerms": ["archive", "container", "package", "storage"]}, {"title": "fas fa-box-open", "searchTerms": ["archive", "container", "package", "storage", "unpack"]}, {"title": "fas fa-boxes", "searchTerms": ["archives", "inventory", "storage", "warehouse"]}, {"title": "fas fa-braille", "searchTerms": ["alphabet", "blind", "dots", "raised", "vision"]}, {"title": "fas fa-brain", "searchTerms": ["cerebellum", "gray matter", "intellect", "medulla oblongata", "mind", "noodle", "wit"]}, {"title": "fas fa-bread-slice", "searchTerms": ["bake", "bakery", "baking", "dough", "flour", "gluten", "grain", "sandwich", "sourdough", "toast", "wheat", "yeast"]}, {"title": "fas fa-briefcase", "searchTerms": ["bag", "business", "luggage", "office", "work"]}, {"title": "fas fa-briefcase-medical", "searchTerms": ["doctor", "emt", "first aid", "health"]}, {"title": "fas fa-broadcast-tower", "searchTerms": ["airwaves", "antenna", "radio", "reception", "waves"]}, {"title": "fas fa-broom", "searchTerms": ["clean", "firebolt", "fly", "halloween", "nimbus 2000", "quidditch", "sweep", "witch"]}, {"title": "fas fa-brush", "searchTerms": ["art", "bristles", "color", "handle", "paint"]}, {"title": "fab fa-btc", "searchTerms": []}, {"title": "fab fa-buffer", "searchTerms": []}, {"title": "fas fa-bug", "searchTerms": ["beetle", "error", "insect", "report"]}, {"title": "fas fa-building", "searchTerms": ["apartment", "business", "city", "company", "office", "work"]}, {"title": "far fa-building", "searchTerms": ["apartment", "business", "city", "company", "office", "work"]}, {"title": "fas fa-bullhorn", "searchTerms": ["announcement", "broadcast", "louder", "megaphone", "share"]}, {"title": "fas fa-bullseye", "searchTerms": ["archery", "goal", "objective", "target"]}, {"title": "fas fa-burn", "searchTerms": ["caliente", "energy", "fire", "flame", "gas", "heat", "hot"]}, {"title": "fab fa-buromobelexperte", "searchTerms": []}, {"title": "fas fa-bus", "searchTerms": ["public transportation", "transportation", "travel", "vehicle"]}, {"title": "fas fa-bus-alt", "searchTerms": ["mta", "public transportation", "transportation", "travel", "vehicle"]}, {"title": "fas fa-business-time", "searchTerms": ["alarm", "briefcase", "business socks", "clock", "flight of the conchords", "reminder", "wednesday"]}, {"title": "fab fa-buy-n-large", "searchTerms": []}, {"title": "fab fa-buysellads", "searchTerms": []}, {"title": "fas fa-calculator", "searchTerms": ["abacus", "addition", "arithmetic", "counting", "math", "multiplication", "subtraction"]}, {"title": "fas fa-calendar", "searchTerms": ["calendar-o", "date", "event", "schedule", "time", "when"]}, {"title": "far fa-calendar", "searchTerms": ["calendar-o", "date", "event", "schedule", "time", "when"]}, {"title": "fas fa-calendar-alt", "searchTerms": ["calendar", "date", "event", "schedule", "time", "when"]}, {"title": "far fa-calendar-alt", "searchTerms": ["calendar", "date", "event", "schedule", "time", "when"]}, {"title": "fas fa-calendar-check", "searchTerms": ["accept", "agree", "appointment", "confirm", "correct", "date", "done", "event", "ok", "schedule", "select", "success", "tick", "time", "todo", "when"]}, {"title": "far fa-calendar-check", "searchTerms": ["accept", "agree", "appointment", "confirm", "correct", "date", "done", "event", "ok", "schedule", "select", "success", "tick", "time", "todo", "when"]}, {"title": "fas fa-calendar-day", "searchTerms": ["date", "detail", "event", "focus", "schedule", "single day", "time", "today", "when"]}, {"title": "fas fa-calendar-minus", "searchTerms": ["calendar", "date", "delete", "event", "negative", "remove", "schedule", "time", "when"]}, {"title": "far fa-calendar-minus", "searchTerms": ["calendar", "date", "delete", "event", "negative", "remove", "schedule", "time", "when"]}, {"title": "fas fa-calendar-plus", "searchTerms": ["add", "calendar", "create", "date", "event", "new", "positive", "schedule", "time", "when"]}, {"title": "far fa-calendar-plus", "searchTerms": ["add", "calendar", "create", "date", "event", "new", "positive", "schedule", "time", "when"]}, {"title": "fas fa-calendar-times", "searchTerms": ["archive", "calendar", "date", "delete", "event", "remove", "schedule", "time", "when", "x"]}, {"title": "far fa-calendar-times", "searchTerms": ["archive", "calendar", "date", "delete", "event", "remove", "schedule", "time", "when", "x"]}, {"title": "fas fa-calendar-week", "searchTerms": ["date", "detail", "event", "focus", "schedule", "single week", "time", "today", "when"]}, {"title": "fas fa-camera", "searchTerms": ["image", "lens", "photo", "picture", "record", "shutter", "video"]}, {"title": "fas fa-camera-retro", "searchTerms": ["image", "lens", "photo", "picture", "record", "shutter", "video"]}, {"title": "fas fa-campground", "searchTerms": ["camping", "fall", "outdoors", "teepee", "tent", "tipi"]}, {"title": "fab fa-canadian-maple-leaf", "searchTerms": ["canada", "flag", "flora", "nature", "plant"]}, {"title": "fas fa-candy-cane", "searchTerms": ["candy", "christmas", "holiday", "mint", "peppermint", "striped", "xmas"]}, {"title": "fas fa-cannabis", "searchTerms": ["bud", "chronic", "drugs", "endica", "endo", "ganja", "marijuana", "mary jane", "pot", "reefer", "sativa", "spliff", "weed", "whacky-tabacky"]}, {"title": "fas fa-capsules", "searchTerms": ["drugs", "medicine", "pills", "prescription"]}, {"title": "fas fa-car", "searchTerms": ["auto", "automobile", "sedan", "transportation", "travel", "vehicle"]}, {"title": "fas fa-car-alt", "searchTerms": ["auto", "automobile", "sedan", "transportation", "travel", "vehicle"]}, {"title": "fas fa-car-battery", "searchTerms": ["auto", "electric", "mechanic", "power"]}, {"title": "fas fa-car-crash", "searchTerms": ["accident", "auto", "automobile", "insurance", "sedan", "transportation", "vehicle", "wreck"]}, {"title": "fas fa-car-side", "searchTerms": ["auto", "automobile", "sedan", "transportation", "travel", "vehicle"]}, {"title": "fas fa-caravan", "searchTerms": ["camper", "motor home", "rv", "trailer", "travel"]}, {"title": "fas fa-caret-down", "searchTerms": ["arrow", "dropdown", "expand", "menu", "more", "triangle"]}, {"title": "fas fa-caret-left", "searchTerms": ["arrow", "back", "previous", "triangle"]}, {"title": "fas fa-caret-right", "searchTerms": ["arrow", "forward", "next", "triangle"]}, {"title": "fas fa-caret-square-down", "searchTerms": ["arrow", "caret-square-o-down", "dropdown", "expand", "menu", "more", "triangle"]}, {"title": "far fa-caret-square-down", "searchTerms": ["arrow", "caret-square-o-down", "dropdown", "expand", "menu", "more", "triangle"]}, {"title": "fas fa-caret-square-left", "searchTerms": ["arrow", "back", "caret-square-o-left", "previous", "triangle"]}, {"title": "far fa-caret-square-left", "searchTerms": ["arrow", "back", "caret-square-o-left", "previous", "triangle"]}, {"title": "fas fa-caret-square-right", "searchTerms": ["arrow", "caret-square-o-right", "forward", "next", "triangle"]}, {"title": "far fa-caret-square-right", "searchTerms": ["arrow", "caret-square-o-right", "forward", "next", "triangle"]}, {"title": "fas fa-caret-square-up", "searchTerms": ["arrow", "caret-square-o-up", "collapse", "triangle", "upload"]}, {"title": "far fa-caret-square-up", "searchTerms": ["arrow", "caret-square-o-up", "collapse", "triangle", "upload"]}, {"title": "fas fa-caret-up", "searchTerms": ["arrow", "collapse", "triangle"]}, {"title": "fas fa-carrot", "searchTerms": ["bugs bunny", "orange", "vegan", "vegetable"]}, {"title": "fas fa-cart-arrow-down", "searchTerms": ["download", "save", "shopping"]}, {"title": "fas fa-cart-plus", "searchTerms": ["add", "create", "new", "positive", "shopping"]}, {"title": "fas fa-cash-register", "searchTerms": ["buy", "cha-ching", "change", "checkout", "commerce", "leaerboard", "machine", "pay", "payment", "purchase", "store"]}, {"title": "fas fa-cat", "searchTerms": ["feline", "halloween", "holiday", "kitten", "kitty", "meow", "pet"]}, {"title": "fab fa-cc-amazon-pay", "searchTerms": []}, {"title": "fab fa-cc-amex", "searchTerms": ["amex"]}, {"title": "fab fa-cc-apple-pay", "searchTerms": []}, {"title": "fab fa-cc-diners-club", "searchTerms": []}, {"title": "fab fa-cc-discover", "searchTerms": []}, {"title": "fab fa-cc-jcb", "searchTerms": []}, {"title": "fab fa-cc-mastercard", "searchTerms": []}, {"title": "fab fa-cc-paypal", "searchTerms": []}, {"title": "fab fa-cc-stripe", "searchTerms": []}, {"title": "fab fa-cc-visa", "searchTerms": []}, {"title": "fab fa-centercode", "searchTerms": []}, {"title": "fab fa-centos", "searchTerms": ["linux", "operating system", "os"]}, {"title": "fas fa-certificate", "searchTerms": ["badge", "star", "verified"]}, {"title": "fas fa-chair", "searchTerms": ["furniture", "seat", "sit"]}, {"title": "fas fa-chalkboard", "searchTerms": ["blackboard", "learning", "school", "teaching", "whiteboard", "writing"]}, {"title": "fas fa-chalkboard-teacher", "searchTerms": ["blackboard", "instructor", "learning", "professor", "school", "whiteboard", "writing"]}, {"title": "fas fa-charging-station", "searchTerms": ["electric", "ev", "tesla", "vehicle"]}, {"title": "fas fa-chart-area", "searchTerms": ["analytics", "area", "chart", "graph"]}, {"title": "fas fa-chart-bar", "searchTerms": ["analytics", "bar", "chart", "graph"]}, {"title": "far fa-chart-bar", "searchTerms": ["analytics", "bar", "chart", "graph"]}, {"title": "fas fa-chart-line", "searchTerms": ["activity", "analytics", "chart", "dashboard", "gain", "graph", "increase", "line"]}, {"title": "fas fa-chart-pie", "searchTerms": ["analytics", "chart", "diagram", "graph", "pie"]}, {"title": "fas fa-check", "searchTerms": ["accept", "agree", "checkmark", "confirm", "correct", "done", "notice", "notification", "notify", "ok", "select", "success", "tick", "todo", "yes"]}, {"title": "fas fa-check-circle", "searchTerms": ["accept", "agree", "confirm", "correct", "done", "ok", "select", "success", "tick", "todo", "yes"]}, {"title": "far fa-check-circle", "searchTerms": ["accept", "agree", "confirm", "correct", "done", "ok", "select", "success", "tick", "todo", "yes"]}, {"title": "fas fa-check-double", "searchTerms": ["accept", "agree", "checkmark", "confirm", "correct", "done", "notice", "notification", "notify", "ok", "select", "success", "tick", "todo"]}, {"title": "fas fa-check-square", "searchTerms": ["accept", "agree", "checkmark", "confirm", "correct", "done", "ok", "select", "success", "tick", "todo", "yes"]}, {"title": "far fa-check-square", "searchTerms": ["accept", "agree", "checkmark", "confirm", "correct", "done", "ok", "select", "success", "tick", "todo", "yes"]}, {"title": "fas fa-cheese", "searchTerms": ["cheddar", "curd", "gouda", "melt", "parmesan", "sandwich", "swiss", "wedge"]}, {"title": "fas fa-chess", "searchTerms": ["board", "castle", "checkmate", "game", "king", "rook", "strategy", "tournament"]}, {"title": "fas fa-chess-bishop", "searchTerms": ["board", "checkmate", "game", "strategy"]}, {"title": "fas fa-chess-board", "searchTerms": ["board", "checkmate", "game", "strategy"]}, {"title": "fas fa-chess-king", "searchTerms": ["board", "checkmate", "game", "strategy"]}, {"title": "fas fa-chess-knight", "searchTerms": ["board", "checkmate", "game", "horse", "strategy"]}, {"title": "fas fa-chess-pawn", "searchTerms": ["board", "checkmate", "game", "strategy"]}, {"title": "fas fa-chess-queen", "searchTerms": ["board", "checkmate", "game", "strategy"]}, {"title": "fas fa-chess-rook", "searchTerms": ["board", "castle", "checkmate", "game", "strategy"]}, {"title": "fas fa-chevron-circle-down", "searchTerms": ["arrow", "download", "dropdown", "menu", "more"]}, {"title": "fas fa-chevron-circle-left", "searchTerms": ["arrow", "back", "previous"]}, {"title": "fas fa-chevron-circle-right", "searchTerms": ["arrow", "forward", "next"]}, {"title": "fas fa-chevron-circle-up", "searchTerms": ["arrow", "collapse", "upload"]}, {"title": "fas fa-chevron-down", "searchTerms": ["arrow", "download", "expand"]}, {"title": "fas fa-chevron-left", "searchTerms": ["arrow", "back", "bracket", "previous"]}, {"title": "fas fa-chevron-right", "searchTerms": ["arrow", "bracket", "forward", "next"]}, {"title": "fas fa-chevron-up", "searchTerms": ["arrow", "collapse", "upload"]}, {"title": "fas fa-child", "searchTerms": ["boy", "girl", "kid", "toddler", "young"]}, {"title": "fab fa-chrome", "searchTerms": ["browser"]}, {"title": "fab fa-chromecast", "searchTerms": []}, {"title": "fas fa-church", "searchTerms": ["building", "cathedral", "chapel", "community", "religion"]}, {"title": "fas fa-circle", "searchTerms": ["circle-thin", "diameter", "dot", "ellipse", "notification", "round"]}, {"title": "far fa-circle", "searchTerms": ["circle-thin", "diameter", "dot", "ellipse", "notification", "round"]}, {"title": "fas fa-circle-notch", "searchTerms": ["circle-o-notch", "diameter", "dot", "ellipse", "round", "spinner"]}, {"title": "fas fa-city", "searchTerms": ["buildings", "busy", "skyscrapers", "urban", "windows"]}, {"title": "fas fa-clinic-medical", "searchTerms": ["doctor", "general practitioner", "hospital", "infirmary", "medicine", "office", "outpatient"]}, {"title": "fas fa-clipboard", "searchTerms": ["copy", "notes", "paste", "record"]}, {"title": "far fa-clipboard", "searchTerms": ["copy", "notes", "paste", "record"]}, {"title": "fas fa-clipboard-check", "searchTerms": ["accept", "agree", "confirm", "done", "ok", "select", "success", "tick", "todo", "yes"]}, {"title": "fas fa-clipboard-list", "searchTerms": ["checklist", "completed", "done", "finished", "intinerary", "ol", "schedule", "tick", "todo", "ul"]}, {"title": "fas fa-clock", "searchTerms": ["date", "late", "schedule", "time", "timer", "timestamp", "watch"]}, {"title": "far fa-clock", "searchTerms": ["date", "late", "schedule", "time", "timer", "timestamp", "watch"]}, {"title": "fas fa-clone", "searchTerms": ["arrange", "copy", "duplicate", "paste"]}, {"title": "far fa-clone", "searchTerms": ["arrange", "copy", "duplicate", "paste"]}, {"title": "fas fa-closed-captioning", "searchTerms": ["cc", "deaf", "hearing", "subtitle", "subtitling", "text", "video"]}, {"title": "far fa-closed-captioning", "searchTerms": ["cc", "deaf", "hearing", "subtitle", "subtitling", "text", "video"]}, {"title": "fas fa-cloud", "searchTerms": ["atmosphere", "fog", "overcast", "save", "upload", "weather"]}, {"title": "fas fa-cloud-download-alt", "searchTerms": ["download", "export", "save"]}, {"title": "fas fa-cloud-meatball", "searchTerms": ["FLDSMDFR", "food", "spaghetti", "storm"]}, {"title": "fas fa-cloud-moon", "searchTerms": ["crescent", "evening", "lunar", "night", "partly cloudy", "sky"]}, {"title": "fas fa-cloud-moon-rain", "searchTerms": ["crescent", "evening", "lunar", "night", "partly cloudy", "precipitation", "rain", "sky", "storm"]}, {"title": "fas fa-cloud-rain", "searchTerms": ["precipitation", "rain", "sky", "storm"]}, {"title": "fas fa-cloud-showers-heavy", "searchTerms": ["precipitation", "rain", "sky", "storm"]}, {"title": "fas fa-cloud-sun", "searchTerms": ["clear", "day", "daytime", "fall", "outdoors", "overcast", "partly cloudy"]}, {"title": "fas fa-cloud-sun-rain", "searchTerms": ["day", "overcast", "precipitation", "storm", "summer", "sunshower"]}, {"title": "fas fa-cloud-upload-alt", "searchTerms": ["cloud-upload", "import", "save", "upload"]}, {"title": "fab fa-cloudscale", "searchTerms": []}, {"title": "fab fa-cloudsmith", "searchTerms": []}, {"title": "fab fa-cloudversify", "searchTerms": []}, {"title": "fas fa-cocktail", "searchTerms": ["alcohol", "beverage", "drink", "gin", "glass", "margarita", "martini", "vodka"]}, {"title": "fas fa-code", "searchTerms": ["brackets", "code", "development", "html"]}, {"title": "fas fa-code-branch", "searchTerms": ["branch", "code-fork", "fork", "git", "github", "rebase", "svn", "vcs", "version"]}, {"title": "fab fa-codepen", "searchTerms": []}, {"title": "fab fa-codiepie", "searchTerms": []}, {"title": "fas fa-coffee", "searchTerms": ["beverage", "breakfast", "cafe", "drink", "fall", "morning", "mug", "seasonal", "tea"]}, {"title": "fas fa-cog", "searchTerms": ["gear", "mechanical", "settings", "sprocket", "wheel"]}, {"title": "fas fa-cogs", "searchTerms": ["gears", "mechanical", "settings", "sprocket", "wheel"]}, {"title": "fas fa-coins", "searchTerms": ["currency", "dime", "financial", "gold", "money", "penny"]}, {"title": "fas fa-columns", "searchTerms": ["browser", "dashboard", "organize", "panes", "split"]}, {"title": "fas fa-comment", "searchTerms": ["bubble", "chat", "commenting", "conversation", "feedback", "message", "note", "notification", "sms", "speech", "texting"]}, {"title": "far fa-comment", "searchTerms": ["bubble", "chat", "commenting", "conversation", "feedback", "message", "note", "notification", "sms", "speech", "texting"]}, {"title": "fas fa-comment-alt", "searchTerms": ["bubble", "chat", "commenting", "conversation", "feedback", "message", "note", "notification", "sms", "speech", "texting"]}, {"title": "far fa-comment-alt", "searchTerms": ["bubble", "chat", "commenting", "conversation", "feedback", "message", "note", "notification", "sms", "speech", "texting"]}, {"title": "fas fa-comment-dollar", "searchTerms": ["bubble", "chat", "commenting", "conversation", "feedback", "message", "money", "note", "notification", "pay", "sms", "speech", "spend", "texting", "transfer"]}, {"title": "fas fa-comment-dots", "searchTerms": ["bubble", "chat", "commenting", "conversation", "feedback", "message", "more", "note", "notification", "reply", "sms", "speech", "texting"]}, {"title": "far fa-comment-dots", "searchTerms": ["bubble", "chat", "commenting", "conversation", "feedback", "message", "more", "note", "notification", "reply", "sms", "speech", "texting"]}, {"title": "fas fa-comment-medical", "searchTerms": ["advice", "bubble", "chat", "commenting", "conversation", "diagnose", "feedback", "message", "note", "notification", "prescription", "sms", "speech", "texting"]}, {"title": "fas fa-comment-slash", "searchTerms": ["bubble", "cancel", "chat", "commenting", "conversation", "feedback", "message", "mute", "note", "notification", "quiet", "sms", "speech", "texting"]}, {"title": "fas fa-comments", "searchTerms": ["bubble", "chat", "commenting", "conversation", "feedback", "message", "note", "notification", "sms", "speech", "texting"]}, {"title": "far fa-comments", "searchTerms": ["bubble", "chat", "commenting", "conversation", "feedback", "message", "note", "notification", "sms", "speech", "texting"]}, {"title": "fas fa-comments-dollar", "searchTerms": ["bubble", "chat", "commenting", "conversation", "feedback", "message", "money", "note", "notification", "pay", "sms", "speech", "spend", "texting", "transfer"]}, {"title": "fas fa-compact-disc", "searchTerms": ["album", "bluray", "cd", "disc", "dvd", "media", "movie", "music", "record", "video", "vinyl"]}, {"title": "fas fa-compass", "searchTerms": ["directions", "directory", "location", "menu", "navigation", "safari", "travel"]}, {"title": "far fa-compass", "searchTerms": ["directions", "directory", "location", "menu", "navigation", "safari", "travel"]}, {"title": "fas fa-compress", "searchTerms": ["collapse", "fullscreen", "minimize", "move", "resize", "shrink", "smaller"]}, {"title": "fas fa-compress-alt", "searchTerms": ["collapse", "fullscreen", "minimize", "move", "resize", "shrink", "smaller"]}, {"title": "fas fa-compress-arrows-alt", "searchTerms": ["collapse", "fullscreen", "minimize", "move", "resize", "shrink", "smaller"]}, {"title": "fas fa-concierge-bell", "searchTerms": ["attention", "hotel", "receptionist", "service", "support"]}, {"title": "fab fa-confluence", "searchTerms": ["atlassian"]}, {"title": "fab fa-connectdevelop", "searchTerms": []}, {"title": "fab fa-contao", "searchTerms": []}, {"title": "fas fa-cookie", "searchTerms": ["baked good", "chips", "chocolate", "eat", "snack", "sweet", "treat"]}, {"title": "fas fa-cookie-bite", "searchTerms": ["baked good", "bitten", "chips", "chocolate", "eat", "snack", "sweet", "treat"]}, {"title": "fas fa-copy", "searchTerms": ["clone", "duplicate", "file", "files-o", "paper", "paste"]}, {"title": "far fa-copy", "searchTerms": ["clone", "duplicate", "file", "files-o", "paper", "paste"]}, {"title": "fas fa-copyright", "searchTerms": ["brand", "mark", "register", "trademark"]}, {"title": "far fa-copyright", "searchTerms": ["brand", "mark", "register", "trademark"]}, {"title": "fab fa-cotton-bureau", "searchTerms": ["clothing", "t-shirts", "tshirts"]}, {"title": "fas fa-couch", "searchTerms": ["chair", "cushion", "furniture", "relax", "sofa"]}, {"title": "fab fa-cpanel", "searchTerms": []}, {"title": "fab fa-creative-commons", "searchTerms": []}, {"title": "fab fa-creative-commons-by", "searchTerms": []}, {"title": "fab fa-creative-commons-nc", "searchTerms": []}, {"title": "fab fa-creative-commons-nc-eu", "searchTerms": []}, {"title": "fab fa-creative-commons-nc-jp", "searchTerms": []}, {"title": "fab fa-creative-commons-nd", "searchTerms": []}, {"title": "fab fa-creative-commons-pd", "searchTerms": []}, {"title": "fab fa-creative-commons-pd-alt", "searchTerms": []}, {"title": "fab fa-creative-commons-remix", "searchTerms": []}, {"title": "fab fa-creative-commons-sa", "searchTerms": []}, {"title": "fab fa-creative-commons-sampling", "searchTerms": []}, {"title": "fab fa-creative-commons-sampling-plus", "searchTerms": []}, {"title": "fab fa-creative-commons-share", "searchTerms": []}, {"title": "fab fa-creative-commons-zero", "searchTerms": []}, {"title": "fas fa-credit-card", "searchTerms": ["buy", "checkout", "credit-card-alt", "debit", "money", "payment", "purchase"]}, {"title": "far fa-credit-card", "searchTerms": ["buy", "checkout", "credit-card-alt", "debit", "money", "payment", "purchase"]}, {"title": "fab fa-critical-role", "searchTerms": ["Dungeons & Dragons", "d&d", "dnd", "fantasy", "game", "gaming", "tabletop"]}, {"title": "fas fa-crop", "searchTerms": ["design", "frame", "mask", "resize", "shrink"]}, {"title": "fas fa-crop-alt", "searchTerms": ["design", "frame", "mask", "resize", "shrink"]}, {"title": "fas fa-cross", "searchTerms": ["catholicism", "christianity", "church", "jesus"]}, {"title": "fas fa-crosshairs", "searchTerms": ["aim", "bullseye", "gpd", "picker", "position"]}, {"title": "fas fa-crow", "searchTerms": ["bird", "bullfrog", "fauna", "halloween", "holiday", "toad"]}, {"title": "fas fa-crown", "searchTerms": ["award", "favorite", "king", "queen", "royal", "tiara"]}, {"title": "fas fa-crutch", "searchTerms": ["cane", "injury", "mobility", "wheelchair"]}, {"title": "fab fa-css3", "searchTerms": ["code"]}, {"title": "fab fa-css3-alt", "searchTerms": []}, {"title": "fas fa-cube", "searchTerms": ["3d", "block", "dice", "package", "square", "tesseract"]}, {"title": "fas fa-cubes", "searchTerms": ["3d", "block", "dice", "package", "pyramid", "square", "stack", "tesseract"]}, {"title": "fas fa-cut", "searchTerms": ["clip", "scissors", "snip"]}, {"title": "fab fa-cuttlefish", "searchTerms": []}, {"title": "fab fa-d-and-d", "searchTerms": []}, {"title": "fab fa-d-and-d-beyond", "searchTerms": ["Dungeons & Dragons", "d&d", "dnd", "fantasy", "gaming", "tabletop"]}, {"title": "fab fa-dailymotion", "searchTerms": []}, {"title": "fab fa-dashcube", "searchTerms": []}, {"title": "fas fa-database", "searchTerms": ["computer", "development", "directory", "memory", "storage"]}, {"title": "fas fa-deaf", "searchTerms": ["ear", "hearing", "sign language"]}, {"title": "fab fa-delicious", "searchTerms": []}, {"title": "fas fa-democrat", "searchTerms": ["american", "democratic party", "donkey", "election", "left", "left-wing", "liberal", "politics", "usa"]}, {"title": "fab fa-deploydog", "searchTerms": []}, {"title": "fab fa-deskpro", "searchTerms": []}, {"title": "fas fa-desktop", "searchTerms": ["computer", "cpu", "demo", "desktop", "device", "imac", "machine", "monitor", "pc", "screen"]}, {"title": "fab fa-dev", "searchTerms": []}, {"title": "fab fa-deviantart", "searchTerms": []}, {"title": "fas fa-dharmachakra", "searchTerms": ["buddhism", "buddhist", "wheel of dharma"]}, {"title": "fab fa-dhl", "searchTerms": ["Dalsey", "Hillblom and Lynn", "german", "package", "shipping"]}, {"title": "fas fa-diagnoses", "searchTerms": ["analyze", "detect", "diagnosis", "examine", "medicine"]}, {"title": "fab fa-diaspora", "searchTerms": []}, {"title": "fas fa-dice", "searchTerms": ["chance", "gambling", "game", "roll"]}, {"title": "fas fa-dice-d20", "searchTerms": ["Dungeons & Dragons", "chance", "d&d", "dnd", "fantasy", "gambling", "game", "roll"]}, {"title": "fas fa-dice-d6", "searchTerms": ["Dungeons & Dragons", "chance", "d&d", "dnd", "fantasy", "gambling", "game", "roll"]}, {"title": "fas fa-dice-five", "searchTerms": ["chance", "gambling", "game", "roll"]}, {"title": "fas fa-dice-four", "searchTerms": ["chance", "gambling", "game", "roll"]}, {"title": "fas fa-dice-one", "searchTerms": ["chance", "gambling", "game", "roll"]}, {"title": "fas fa-dice-six", "searchTerms": ["chance", "gambling", "game", "roll"]}, {"title": "fas fa-dice-three", "searchTerms": ["chance", "gambling", "game", "roll"]}, {"title": "fas fa-dice-two", "searchTerms": ["chance", "gambling", "game", "roll"]}, {"title": "fab fa-digg", "searchTerms": []}, {"title": "fab fa-digital-ocean", "searchTerms": []}, {"title": "fas fa-digital-tachograph", "searchTerms": ["data", "distance", "speed", "tachometer"]}, {"title": "fas fa-directions", "searchTerms": ["map", "navigation", "sign", "turn"]}, {"title": "fab fa-discord", "searchTerms": []}, {"title": "fab fa-discourse", "searchTerms": []}, {"title": "fas fa-divide", "searchTerms": ["arithmetic", "calculus", "division", "math"]}, {"title": "fas fa-dizzy", "searchTerms": ["dazed", "dead", "disapprove", "emoticon", "face"]}, {"title": "far fa-dizzy", "searchTerms": ["dazed", "dead", "disapprove", "emoticon", "face"]}, {"title": "fas fa-dna", "searchTerms": ["double helix", "genetic", "helix", "molecule", "protein"]}, {"title": "fab fa-dochub", "searchTerms": []}, {"title": "fab fa-docker", "searchTerms": []}, {"title": "fas fa-dog", "searchTerms": ["animal", "canine", "fauna", "mammal", "pet", "pooch", "puppy", "woof"]}, {"title": "fas fa-dollar-sign", "searchTerms": ["$", "cost", "dollar-sign", "money", "price", "usd"]}, {"title": "fas fa-dolly", "searchTerms": ["carry", "shipping", "transport"]}, {"title": "fas fa-dolly-flatbed", "searchTerms": ["carry", "inventory", "shipping", "transport"]}, {"title": "fas fa-donate", "searchTerms": ["contribute", "generosity", "gift", "give"]}, {"title": "fas fa-door-closed", "searchTerms": ["enter", "exit", "locked"]}, {"title": "fas fa-door-open", "searchTerms": ["enter", "exit", "welcome"]}, {"title": "fas fa-dot-circle", "searchTerms": ["bullseye", "notification", "target"]}, {"title": "far fa-dot-circle", "searchTerms": ["bullseye", "notification", "target"]}, {"title": "fas fa-dove", "searchTerms": ["bird", "fauna", "flying", "peace", "war"]}, {"title": "fas fa-download", "searchTerms": ["export", "hard drive", "save", "transfer"]}, {"title": "fab fa-draft2digital", "searchTerms": []}, {"title": "fas fa-drafting-compass", "searchTerms": ["design", "map", "mechanical drawing", "plot", "plotting"]}, {"title": "fas fa-dragon", "searchTerms": ["Dungeons & Dragons", "d&d", "dnd", "fantasy", "fire", "lizard", "serpent"]}, {"title": "fas fa-draw-polygon", "searchTerms": ["anchors", "lines", "object", "render", "shape"]}, {"title": "fab fa-dribbble", "searchTerms": []}, {"title": "fab fa-dribbble-square", "searchTerms": []}, {"title": "fab fa-dropbox", "searchTerms": []}, {"title": "fas fa-drum", "searchTerms": ["instrument", "music", "percussion", "snare", "sound"]}, {"title": "fas fa-drum-steelpan", "searchTerms": ["calypso", "instrument", "music", "percussion", "reggae", "snare", "sound", "steel", "tropical"]}, {"title": "fas fa-drumstick-bite", "searchTerms": ["bone", "chicken", "leg", "meat", "poultry", "turkey"]}, {"title": "fab fa-drupal", "searchTerms": []}, {"title": "fas fa-dumbbell", "searchTerms": ["exercise", "gym", "strength", "weight", "weight-lifting"]}, {"title": "fas fa-dumpster", "searchTerms": ["alley", "bin", "commercial", "trash", "waste"]}, {"title": "fas fa-dumpster-fire", "searchTerms": ["alley", "bin", "commercial", "danger", "dangerous", "euphemism", "flame", "heat", "hot", "trash", "waste"]}, {"title": "fas fa-dungeon", "searchTerms": ["Dungeons & Dragons", "building", "d&d", "dnd", "door", "entrance", "fantasy", "gate"]}, {"title": "fab fa-dyalog", "searchTerms": []}, {"title": "fab fa-earlybirds", "searchTerms": []}, {"title": "fab fa-ebay", "searchTerms": []}, {"title": "fab fa-edge", "searchTerms": ["browser", "ie"]}, {"title": "fas fa-edit", "searchTerms": ["edit", "pen", "pencil", "update", "write"]}, {"title": "far fa-edit", "searchTerms": ["edit", "pen", "pencil", "update", "write"]}, {"title": "fas fa-egg", "searchTerms": ["breakfast", "chicken", "easter", "shell", "yolk"]}, {"title": "fas fa-eject", "searchTerms": ["abort", "cancel", "cd", "discharge"]}, {"title": "fab fa-elementor", "searchTerms": []}, {"title": "fas fa-ellipsis-h", "searchTerms": ["dots", "drag", "kebab", "list", "menu", "nav", "navigation", "ol", "reorder", "settings", "ul"]}, {"title": "fas fa-ellipsis-v", "searchTerms": ["dots", "drag", "kebab", "list", "menu", "nav", "navigation", "ol", "reorder", "settings", "ul"]}, {"title": "fab fa-ello", "searchTerms": []}, {"title": "fab fa-ember", "searchTerms": []}, {"title": "fab fa-empire", "searchTerms": []}, {"title": "fas fa-envelope", "searchTerms": ["e-mail", "email", "letter", "mail", "message", "notification", "support"]}, {"title": "far fa-envelope", "searchTerms": ["e-mail", "email", "letter", "mail", "message", "notification", "support"]}, {"title": "fas fa-envelope-open", "searchTerms": ["e-mail", "email", "letter", "mail", "message", "notification", "support"]}, {"title": "far fa-envelope-open", "searchTerms": ["e-mail", "email", "letter", "mail", "message", "notification", "support"]}, {"title": "fas fa-envelope-open-text", "searchTerms": ["e-mail", "email", "letter", "mail", "message", "notification", "support"]}, {"title": "fas fa-envelope-square", "searchTerms": ["e-mail", "email", "letter", "mail", "message", "notification", "support"]}, {"title": "fab fa-envira", "searchTerms": ["leaf"]}, {"title": "fas fa-equals", "searchTerms": ["arithmetic", "even", "match", "math"]}, {"title": "fas fa-eraser", "searchTerms": ["art", "delete", "remove", "rubber"]}, {"title": "fab fa-erlang", "searchTerms": []}, {"title": "fab fa-ethereum", "searchTerms": []}, {"title": "fas fa-ethernet", "searchTerms": ["cable", "cat 5", "cat 6", "connection", "hardware", "internet", "network", "wired"]}, {"title": "fab fa-etsy", "searchTerms": []}, {"title": "fas fa-euro-sign", "searchTerms": ["currency", "dollar", "exchange", "money"]}, {"title": "fab fa-evernote", "searchTerms": []}, {"title": "fas fa-exchange-alt", "searchTerms": ["arrow", "arrows", "exchange", "reciprocate", "return", "swap", "transfer"]}, {"title": "fas fa-exclamation", "searchTerms": ["alert", "danger", "error", "important", "notice", "notification", "notify", "problem", "warning"]}, {"title": "fas fa-exclamation-circle", "searchTerms": ["alert", "danger", "error", "important", "notice", "notification", "notify", "problem", "warning"]}, {"title": "fas fa-exclamation-triangle", "searchTerms": ["alert", "danger", "error", "important", "notice", "notification", "notify", "problem", "warning"]}, {"title": "fas fa-expand", "searchTerms": ["arrow", "bigger", "enlarge", "resize"]}, {"title": "fas fa-expand-alt", "searchTerms": ["arrow", "bigger", "enlarge", "resize"]}, {"title": "fas fa-expand-arrows-alt", "searchTerms": ["arrows-alt", "bigger", "enlarge", "move", "resize"]}, {"title": "fab fa-expeditedssl", "searchTerms": []}, {"title": "fas fa-external-link-alt", "searchTerms": ["external-link", "new", "open", "share"]}, {"title": "fas fa-external-link-square-alt", "searchTerms": ["external-link-square", "new", "open", "share"]}, {"title": "fas fa-eye", "searchTerms": ["look", "optic", "see", "seen", "show", "sight", "views", "visible"]}, {"title": "far fa-eye", "searchTerms": ["look", "optic", "see", "seen", "show", "sight", "views", "visible"]}, {"title": "fas fa-eye-dropper", "searchTerms": ["beaker", "clone", "color", "copy", "eyedropper", "pipette"]}, {"title": "fas fa-eye-slash", "searchTerms": ["blind", "hide", "show", "toggle", "unseen", "views", "visible", "visiblity"]}, {"title": "far fa-eye-slash", "searchTerms": ["blind", "hide", "show", "toggle", "unseen", "views", "visible", "visiblity"]}, {"title": "fab fa-facebook", "searchTerms": ["facebook-official", "social network"]}, {"title": "fab fa-facebook-f", "searchTerms": ["facebook"]}, {"title": "fab fa-facebook-messenger", "searchTerms": []}, {"title": "fab fa-facebook-square", "searchTerms": ["social network"]}, {"title": "fas fa-fan", "searchTerms": ["ac", "air conditioning", "blade", "blower", "cool", "hot"]}, {"title": "fab fa-fantasy-flight-games", "searchTerms": ["Dungeons & Dragons", "d&d", "dnd", "fantasy", "game", "gaming", "tabletop"]}, {"title": "fas fa-fast-backward", "searchTerms": ["beginning", "first", "previous", "rewind", "start"]}, {"title": "fas fa-fast-forward", "searchTerms": ["end", "last", "next"]}, {"title": "fas fa-fax", "searchTerms": ["business", "communicate", "copy", "facsimile", "send"]}, {"title": "fas fa-feather", "searchTerms": ["bird", "light", "plucked", "quill", "write"]}, {"title": "fas fa-feather-alt", "searchTerms": ["bird", "light", "plucked", "quill", "write"]}, {"title": "fab fa-fedex", "searchTerms": ["Federal Express", "package", "shipping"]}, {"title": "fab fa-fedora", "searchTerms": ["linux", "operating system", "os"]}, {"title": "fas fa-female", "searchTerms": ["human", "person", "profile", "user", "woman"]}, {"title": "fas fa-fighter-jet", "searchTerms": ["airplane", "fast", "fly", "goose", "maverick", "plane", "quick", "top gun", "transportation", "travel"]}, {"title": "fab fa-figma", "searchTerms": ["app", "design", "interface"]}, {"title": "fas fa-file", "searchTerms": ["document", "new", "page", "pdf", "resume"]}, {"title": "far fa-file", "searchTerms": ["document", "new", "page", "pdf", "resume"]}, {"title": "fas fa-file-alt", "searchTerms": ["document", "file-text", "invoice", "new", "page", "pdf"]}, {"title": "far fa-file-alt", "searchTerms": ["document", "file-text", "invoice", "new", "page", "pdf"]}, {"title": "fas fa-file-archive", "searchTerms": [".zip", "bundle", "compress", "compression", "download", "zip"]}, {"title": "far fa-file-archive", "searchTerms": [".zip", "bundle", "compress", "compression", "download", "zip"]}, {"title": "fas fa-file-audio", "searchTerms": ["document", "mp3", "music", "page", "play", "sound"]}, {"title": "far fa-file-audio", "searchTerms": ["document", "mp3", "music", "page", "play", "sound"]}, {"title": "fas fa-file-code", "searchTerms": ["css", "development", "document", "html"]}, {"title": "far fa-file-code", "searchTerms": ["css", "development", "document", "html"]}, {"title": "fas fa-file-contract", "searchTerms": ["agreement", "binding", "document", "legal", "signature"]}, {"title": "fas fa-file-csv", "searchTerms": ["document", "excel", "numbers", "spreadsheets", "table"]}, {"title": "fas fa-file-download", "searchTerms": ["document", "export", "save"]}, {"title": "fas fa-file-excel", "searchTerms": ["csv", "document", "numbers", "spreadsheets", "table"]}, {"title": "far fa-file-excel", "searchTerms": ["csv", "document", "numbers", "spreadsheets", "table"]}, {"title": "fas fa-file-export", "searchTerms": ["download", "save"]}, {"title": "fas fa-file-image", "searchTerms": ["document", "image", "jpg", "photo", "png"]}, {"title": "far fa-file-image", "searchTerms": ["document", "image", "jpg", "photo", "png"]}, {"title": "fas fa-file-import", "searchTerms": ["copy", "document", "send", "upload"]}, {"title": "fas fa-file-invoice", "searchTerms": ["account", "bill", "charge", "document", "payment", "receipt"]}, {"title": "fas fa-file-invoice-dollar", "searchTerms": ["$", "account", "bill", "charge", "document", "dollar-sign", "money", "payment", "receipt", "usd"]}, {"title": "fas fa-file-medical", "searchTerms": ["document", "health", "history", "prescription", "record"]}, {"title": "fas fa-file-medical-alt", "searchTerms": ["document", "health", "history", "prescription", "record"]}, {"title": "fas fa-file-pdf", "searchTerms": ["acrobat", "document", "preview", "save"]}, {"title": "far fa-file-pdf", "searchTerms": ["acrobat", "document", "preview", "save"]}, {"title": "fas fa-file-powerpoint", "searchTerms": ["display", "document", "keynote", "presentation"]}, {"title": "far fa-file-powerpoint", "searchTerms": ["display", "document", "keynote", "presentation"]}, {"title": "fas fa-file-prescription", "searchTerms": ["document", "drugs", "medical", "medicine", "rx"]}, {"title": "fas fa-file-signature", "searchTerms": ["John Hancock", "contract", "document", "name"]}, {"title": "fas fa-file-upload", "searchTerms": ["document", "import", "page", "save"]}, {"title": "fas fa-file-video", "searchTerms": ["document", "m4v", "movie", "mp4", "play"]}, {"title": "far fa-file-video", "searchTerms": ["document", "m4v", "movie", "mp4", "play"]}, {"title": "fas fa-file-word", "searchTerms": ["document", "edit", "page", "text", "writing"]}, {"title": "far fa-file-word", "searchTerms": ["document", "edit", "page", "text", "writing"]}, {"title": "fas fa-fill", "searchTerms": ["bucket", "color", "paint", "paint bucket"]}, {"title": "fas fa-fill-drip", "searchTerms": ["bucket", "color", "drop", "paint", "paint bucket", "spill"]}, {"title": "fas fa-film", "searchTerms": ["cinema", "movie", "strip", "video"]}, {"title": "fas fa-filter", "searchTerms": ["funnel", "options", "separate", "sort"]}, {"title": "fas fa-fingerprint", "searchTerms": ["human", "id", "identification", "lock", "smudge", "touch", "unique", "unlock"]}, {"title": "fas fa-fire", "searchTerms": ["burn", "caliente", "flame", "heat", "hot", "popular"]}, {"title": "fas fa-fire-alt", "searchTerms": ["burn", "caliente", "flame", "heat", "hot", "popular"]}, {"title": "fas fa-fire-extinguisher", "searchTerms": ["burn", "caliente", "fire fighter", "flame", "heat", "hot", "rescue"]}, {"title": "fab fa-firefox", "searchTerms": ["browser"]}, {"title": "fab fa-firefox-browser", "searchTerms": ["browser"]}, {"title": "fas fa-first-aid", "searchTerms": ["emergency", "emt", "health", "medical", "rescue"]}, {"title": "fab fa-first-order", "searchTerms": []}, {"title": "fab fa-first-order-alt", "searchTerms": []}, {"title": "fab fa-firstdraft", "searchTerms": []}, {"title": "fas fa-fish", "searchTerms": ["fauna", "gold", "seafood", "swimming"]}, {"title": "fas fa-fist-raised", "searchTerms": ["Dungeons & Dragons", "d&d", "dnd", "fantasy", "hand", "ki", "monk", "resist", "strength", "unarmed combat"]}, {"title": "fas fa-flag", "searchTerms": ["country", "notice", "notification", "notify", "pole", "report", "symbol"]}, {"title": "far fa-flag", "searchTerms": ["country", "notice", "notification", "notify", "pole", "report", "symbol"]}, {"title": "fas fa-flag-checkered", "searchTerms": ["notice", "notification", "notify", "pole", "racing", "report", "symbol"]}, {"title": "fas fa-flag-usa", "searchTerms": ["betsy ross", "country", "old glory", "stars", "stripes", "symbol"]}, {"title": "fas fa-flask", "searchTerms": ["beaker", "experimental", "labs", "science"]}, {"title": "fab fa-flickr", "searchTerms": []}, {"title": "fab fa-flipboard", "searchTerms": []}, {"title": "fas fa-flushed", "searchTerms": ["embarrassed", "emoticon", "face"]}, {"title": "far fa-flushed", "searchTerms": ["embarrassed", "emoticon", "face"]}, {"title": "fab fa-fly", "searchTerms": []}, {"title": "fas fa-folder", "searchTerms": ["archive", "directory", "document", "file"]}, {"title": "far fa-folder", "searchTerms": ["archive", "directory", "document", "file"]}, {"title": "fas fa-folder-minus", "searchTerms": ["archive", "delete", "directory", "document", "file", "negative", "remove"]}, {"title": "fas fa-folder-open", "searchTerms": ["archive", "directory", "document", "empty", "file", "new"]}, {"title": "far fa-folder-open", "searchTerms": ["archive", "directory", "document", "empty", "file", "new"]}, {"title": "fas fa-folder-plus", "searchTerms": ["add", "archive", "create", "directory", "document", "file", "new", "positive"]}, {"title": "fas fa-font", "searchTerms": ["alphabet", "glyph", "text", "type", "typeface"]}, {"title": "fab fa-font-awesome", "searchTerms": ["meanpath"]}, {"title": "fab fa-font-awesome-alt", "searchTerms": []}, {"title": "fab fa-font-awesome-flag", "searchTerms": []}, {"title": "far fa-font-awesome-logo-full", "searchTerms": []}, {"title": "fas fa-font-awesome-logo-full", "searchTerms": []}, {"title": "fab fa-font-awesome-logo-full", "searchTerms": []}, {"title": "fab fa-fonticons", "searchTerms": []}, {"title": "fab fa-fonticons-fi", "searchTerms": []}, {"title": "fas fa-football-ball", "searchTerms": ["ball", "fall", "nfl", "pigskin", "seasonal"]}, {"title": "fab fa-fort-awesome", "searchTerms": ["castle"]}, {"title": "fab fa-fort-awesome-alt", "searchTerms": ["castle"]}, {"title": "fab fa-forumbee", "searchTerms": []}, {"title": "fas fa-forward", "searchTerms": ["forward", "next", "skip"]}, {"title": "fab fa-foursquare", "searchTerms": []}, {"title": "fab fa-free-code-camp", "searchTerms": []}, {"title": "fab fa-freebsd", "searchTerms": []}, {"title": "fas fa-frog", "searchTerms": ["amphibian", "bullfrog", "fauna", "hop", "kermit", "kiss", "prince", "ribbit", "toad", "wart"]}, {"title": "fas fa-frown", "searchTerms": ["disapprove", "emoticon", "face", "rating", "sad"]}, {"title": "far fa-frown", "searchTerms": ["disapprove", "emoticon", "face", "rating", "sad"]}, {"title": "fas fa-frown-open", "searchTerms": ["disapprove", "emoticon", "face", "rating", "sad"]}, {"title": "far fa-frown-open", "searchTerms": ["disapprove", "emoticon", "face", "rating", "sad"]}, {"title": "fab fa-fulcrum", "searchTerms": []}, {"title": "fas fa-funnel-dollar", "searchTerms": ["filter", "money", "options", "separate", "sort"]}, {"title": "fas fa-futbol", "searchTerms": ["ball", "football", "mls", "soccer"]}, {"title": "far fa-futbol", "searchTerms": ["ball", "football", "mls", "soccer"]}, {"title": "fab fa-galactic-republic", "searchTerms": ["politics", "star wars"]}, {"title": "fab fa-galactic-senate", "searchTerms": ["star wars"]}, {"title": "fas fa-gamepad", "searchTerms": ["arcade", "controller", "d-pad", "joystick", "video", "video game"]}, {"title": "fas fa-gas-pump", "searchTerms": ["car", "fuel", "gasoline", "petrol"]}, {"title": "fas fa-gavel", "searchTerms": ["hammer", "judge", "law", "lawyer", "opinion"]}, {"title": "fas fa-gem", "searchTerms": ["diamond", "jewelry", "sapphire", "stone", "treasure"]}, {"title": "far fa-gem", "searchTerms": ["diamond", "jewelry", "sapphire", "stone", "treasure"]}, {"title": "fas fa-genderless", "searchTerms": ["androgynous", "asexual", "sexless"]}, {"title": "fab fa-get-pocket", "searchTerms": []}, {"title": "fab fa-gg", "searchTerms": []}, {"title": "fab fa-gg-circle", "searchTerms": []}, {"title": "fas fa-ghost", "searchTerms": ["apparition", "blinky", "clyde", "floating", "halloween", "holiday", "inky", "pinky", "spirit"]}, {"title": "fas fa-gift", "searchTerms": ["christmas", "generosity", "giving", "holiday", "party", "present", "wrapped", "xmas"]}, {"title": "fas fa-gifts", "searchTerms": ["christmas", "generosity", "giving", "holiday", "party", "present", "wrapped", "xmas"]}, {"title": "fab fa-git", "searchTerms": []}, {"title": "fab fa-git-alt", "searchTerms": []}, {"title": "fab fa-git-square", "searchTerms": []}, {"title": "fab fa-github", "searchTerms": ["octocat"]}, {"title": "fab fa-github-alt", "searchTerms": ["octocat"]}, {"title": "fab fa-github-square", "searchTerms": ["octocat"]}, {"title": "fab fa-gitkraken", "searchTerms": []}, {"title": "fab fa-gitlab", "searchTerms": ["Axosoft"]}, {"title": "fab fa-gitter", "searchTerms": []}, {"title": "fas fa-glass-cheers", "searchTerms": ["alcohol", "bar", "beverage", "celebration", "champagne", "clink", "drink", "holiday", "new year's eve", "party", "toast"]}, {"title": "fas fa-glass-martini", "searchTerms": ["alcohol", "bar", "beverage", "drink", "liquor"]}, {"title": "fas fa-glass-martini-alt", "searchTerms": ["alcohol", "bar", "beverage", "drink", "liquor"]}, {"title": "fas fa-glass-whiskey", "searchTerms": ["alcohol", "bar", "beverage", "bourbon", "drink", "liquor", "neat", "rye", "scotch", "whisky"]}, {"title": "fas fa-glasses", "searchTerms": ["hipster", "nerd", "reading", "sight", "spectacles", "vision"]}, {"title": "fab fa-glide", "searchTerms": []}, {"title": "fab fa-glide-g", "searchTerms": []}, {"title": "fas fa-globe", "searchTerms": ["all", "coordinates", "country", "earth", "global", "gps", "language", "localize", "location", "map", "online", "place", "planet", "translate", "travel", "world"]}, {"title": "fas fa-globe-africa", "searchTerms": ["all", "country", "earth", "global", "gps", "language", "localize", "location", "map", "online", "place", "planet", "translate", "travel", "world"]}, {"title": "fas fa-globe-americas", "searchTerms": ["all", "country", "earth", "global", "gps", "language", "localize", "location", "map", "online", "place", "planet", "translate", "travel", "world"]}, {"title": "fas fa-globe-asia", "searchTerms": ["all", "country", "earth", "global", "gps", "language", "localize", "location", "map", "online", "place", "planet", "translate", "travel", "world"]}, {"title": "fas fa-globe-europe", "searchTerms": ["all", "country", "earth", "global", "gps", "language", "localize", "location", "map", "online", "place", "planet", "translate", "travel", "world"]}, {"title": "fab fa-gofore", "searchTerms": []}, {"title": "fas fa-golf-ball", "searchTerms": ["caddy", "eagle", "putt", "tee"]}, {"title": "fab fa-goodreads", "searchTerms": []}, {"title": "fab fa-goodreads-g", "searchTerms": []}, {"title": "fab fa-google", "searchTerms": []}, {"title": "fab fa-google-drive", "searchTerms": []}, {"title": "fab fa-google-play", "searchTerms": []}, {"title": "fab fa-google-plus", "searchTerms": ["google-plus-circle", "google-plus-official"]}, {"title": "fab fa-google-plus-g", "searchTerms": ["google-plus", "social network"]}, {"title": "fab fa-google-plus-square", "searchTerms": ["social network"]}, {"title": "fab fa-google-wallet", "searchTerms": []}, {"title": "fas fa-gopuram", "searchTerms": ["building", "entrance", "hinduism", "temple", "tower"]}, {"title": "fas fa-graduation-cap", "searchTerms": ["ceremony", "college", "graduate", "learning", "school", "student"]}, {"title": "fab fa-gratipay", "searchTerms": ["favorite", "heart", "like", "love"]}, {"title": "fab fa-grav", "searchTerms": []}, {"title": "fas fa-greater-than", "searchTerms": ["arithmetic", "compare", "math"]}, {"title": "fas fa-greater-than-equal", "searchTerms": ["arithmetic", "compare", "math"]}, {"title": "fas fa-grimace", "searchTerms": ["cringe", "emoticon", "face", "teeth"]}, {"title": "far fa-grimace", "searchTerms": ["cringe", "emoticon", "face", "teeth"]}, {"title": "fas fa-grin", "searchTerms": ["emoticon", "face", "laugh", "smile"]}, {"title": "far fa-grin", "searchTerms": ["emoticon", "face", "laugh", "smile"]}, {"title": "fas fa-grin-alt", "searchTerms": ["emoticon", "face", "laugh", "smile"]}, {"title": "far fa-grin-alt", "searchTerms": ["emoticon", "face", "laugh", "smile"]}, {"title": "fas fa-grin-beam", "searchTerms": ["emoticon", "face", "laugh", "smile"]}, {"title": "far fa-grin-beam", "searchTerms": ["emoticon", "face", "laugh", "smile"]}, {"title": "fas fa-grin-beam-sweat", "searchTerms": ["embarass", "emoticon", "face", "smile"]}, {"title": "far fa-grin-beam-sweat", "searchTerms": ["embarass", "emoticon", "face", "smile"]}, {"title": "fas fa-grin-hearts", "searchTerms": ["emoticon", "face", "love", "smile"]}, {"title": "far fa-grin-hearts", "searchTerms": ["emoticon", "face", "love", "smile"]}, {"title": "fas fa-grin-squint", "searchTerms": ["emoticon", "face", "laugh", "smile"]}, {"title": "far fa-grin-squint", "searchTerms": ["emoticon", "face", "laugh", "smile"]}, {"title": "fas fa-grin-squint-tears", "searchTerms": ["emoticon", "face", "happy", "smile"]}, {"title": "far fa-grin-squint-tears", "searchTerms": ["emoticon", "face", "happy", "smile"]}, {"title": "fas fa-grin-stars", "searchTerms": ["emoticon", "face", "star-struck"]}, {"title": "far fa-grin-stars", "searchTerms": ["emoticon", "face", "star-struck"]}, {"title": "fas fa-grin-tears", "searchTerms": ["LOL", "emoticon", "face"]}, {"title": "far fa-grin-tears", "searchTerms": ["LOL", "emoticon", "face"]}, {"title": "fas fa-grin-tongue", "searchTerms": ["LOL", "emoticon", "face"]}, {"title": "far fa-grin-tongue", "searchTerms": ["LOL", "emoticon", "face"]}, {"title": "fas fa-grin-tongue-squint", "searchTerms": ["LOL", "emoticon", "face"]}, {"title": "far fa-grin-tongue-squint", "searchTerms": ["LOL", "emoticon", "face"]}, {"title": "fas fa-grin-tongue-wink", "searchTerms": ["LOL", "emoticon", "face"]}, {"title": "far fa-grin-tongue-wink", "searchTerms": ["LOL", "emoticon", "face"]}, {"title": "fas fa-grin-wink", "searchTerms": ["emoticon", "face", "flirt", "laugh", "smile"]}, {"title": "far fa-grin-wink", "searchTerms": ["emoticon", "face", "flirt", "laugh", "smile"]}, {"title": "fas fa-grip-horizontal", "searchTerms": ["affordance", "drag", "drop", "grab", "handle"]}, {"title": "fas fa-grip-lines", "searchTerms": ["affordance", "drag", "drop", "grab", "handle"]}, {"title": "fas fa-grip-lines-vertical", "searchTerms": ["affordance", "drag", "drop", "grab", "handle"]}, {"title": "fas fa-grip-vertical", "searchTerms": ["affordance", "drag", "drop", "grab", "handle"]}, {"title": "fab fa-gripfire", "searchTerms": []}, {"title": "fab fa-grunt", "searchTerms": []}, {"title": "fas fa-guitar", "searchTerms": ["acoustic", "instrument", "music", "rock", "rock and roll", "song", "strings"]}, {"title": "fab fa-gulp", "searchTerms": []}, {"title": "fas fa-h-square", "searchTerms": ["directions", "emergency", "hospital", "hotel", "map"]}, {"title": "fab fa-hacker-news", "searchTerms": []}, {"title": "fab fa-hacker-news-square", "searchTerms": []}, {"title": "fab fa-hackerrank", "searchTerms": []}, {"title": "fas fa-hamburger", "searchTerms": ["bacon", "beef", "burger", "burger king", "cheeseburger", "fast food", "grill", "ground beef", "mcdonalds", "sandwich"]}, {"title": "fas fa-hammer", "searchTerms": ["admin", "fix", "repair", "settings", "tool"]}, {"title": "fas fa-hamsa", "searchTerms": ["amulet", "christianity", "islam", "jewish", "judaism", "muslim", "protection"]}, {"title": "fas fa-hand-holding", "searchTerms": ["carry", "lift"]}, {"title": "fas fa-hand-holding-heart", "searchTerms": ["carry", "charity", "gift", "lift", "package"]}, {"title": "fas fa-hand-holding-usd", "searchTerms": ["$", "carry", "dollar sign", "donation", "giving", "lift", "money", "price"]}, {"title": "fas fa-hand-lizard", "searchTerms": ["game", "roshambo"]}, {"title": "far fa-hand-lizard", "searchTerms": ["game", "roshambo"]}, {"title": "fas fa-hand-middle-finger", "searchTerms": ["flip the bird", "gesture", "hate", "rude"]}, {"title": "fas fa-hand-paper", "searchTerms": ["game", "halt", "roshambo", "stop"]}, {"title": "far fa-hand-paper", "searchTerms": ["game", "halt", "roshambo", "stop"]}, {"title": "fas fa-hand-peace", "searchTerms": ["rest", "truce"]}, {"title": "far fa-hand-peace", "searchTerms": ["rest", "truce"]}, {"title": "fas fa-hand-point-down", "searchTerms": ["finger", "hand-o-down", "point"]}, {"title": "far fa-hand-point-down", "searchTerms": ["finger", "hand-o-down", "point"]}, {"title": "fas fa-hand-point-left", "searchTerms": ["back", "finger", "hand-o-left", "left", "point", "previous"]}, {"title": "far fa-hand-point-left", "searchTerms": ["back", "finger", "hand-o-left", "left", "point", "previous"]}, {"title": "fas fa-hand-point-right", "searchTerms": ["finger", "forward", "hand-o-right", "next", "point", "right"]}, {"title": "far fa-hand-point-right", "searchTerms": ["finger", "forward", "hand-o-right", "next", "point", "right"]}, {"title": "fas fa-hand-point-up", "searchTerms": ["finger", "hand-o-up", "point"]}, {"title": "far fa-hand-point-up", "searchTerms": ["finger", "hand-o-up", "point"]}, {"title": "fas fa-hand-pointer", "searchTerms": ["arrow", "cursor", "select"]}, {"title": "far fa-hand-pointer", "searchTerms": ["arrow", "cursor", "select"]}, {"title": "fas fa-hand-rock", "searchTerms": ["fist", "game", "roshambo"]}, {"title": "far fa-hand-rock", "searchTerms": ["fist", "game", "roshambo"]}, {"title": "fas fa-hand-scissors", "searchTerms": ["cut", "game", "roshambo"]}, {"title": "far fa-hand-scissors", "searchTerms": ["cut", "game", "roshambo"]}, {"title": "fas fa-hand-spock", "searchTerms": ["live long", "prosper", "salute", "star trek", "vulcan"]}, {"title": "far fa-hand-spock", "searchTerms": ["live long", "prosper", "salute", "star trek", "vulcan"]}, {"title": "fas fa-hands", "searchTerms": ["carry", "hold", "lift"]}, {"title": "fas fa-hands-helping", "searchTerms": ["aid", "assistance", "handshake", "partnership", "volunteering"]}, {"title": "fas fa-handshake", "searchTerms": ["agreement", "greeting", "meeting", "partnership"]}, {"title": "far fa-handshake", "searchTerms": ["agreement", "greeting", "meeting", "partnership"]}, {"title": "fas fa-hanukiah", "searchTerms": ["candle", "hanukkah", "jewish", "judaism", "light"]}, {"title": "fas fa-hard-hat", "searchTerms": ["construction", "hardhat", "helmet", "safety"]}, {"title": "fas fa-hashtag", "searchTerms": ["Twitter", "instagram", "pound", "social media", "tag"]}, {"title": "fas fa-hat-cowboy", "searchTerms": ["buckaroo", "horse", "jackeroo", "john b.", "old west", "pardner", "ranch", "rancher", "rodeo", "western", "wrangler"]}, {"title": "fas fa-hat-cowboy-side", "searchTerms": ["buckaroo", "horse", "jackeroo", "john b.", "old west", "pardner", "ranch", "rancher", "rodeo", "western", "wrangler"]}, {"title": "fas fa-hat-wizard", "searchTerms": ["Dungeons & Dragons", "accessory", "buckle", "clothing", "d&d", "dnd", "fantasy", "halloween", "head", "holiday", "mage", "magic", "pointy", "witch"]}, {"title": "fas fa-hdd", "searchTerms": ["cpu", "hard drive", "harddrive", "machine", "save", "storage"]}, {"title": "far fa-hdd", "searchTerms": ["cpu", "hard drive", "harddrive", "machine", "save", "storage"]}, {"title": "fas fa-heading", "searchTerms": ["format", "header", "text", "title"]}, {"title": "fas fa-headphones", "searchTerms": ["audio", "listen", "music", "sound", "speaker"]}, {"title": "fas fa-headphones-alt", "searchTerms": ["audio", "listen", "music", "sound", "speaker"]}, {"title": "fas fa-headset", "searchTerms": ["audio", "gamer", "gaming", "listen", "live chat", "microphone", "shot caller", "sound", "support", "telemarketer"]}, {"title": "fas fa-heart", "searchTerms": ["favorite", "like", "love", "relationship", "valentine"]}, {"title": "far fa-heart", "searchTerms": ["favorite", "like", "love", "relationship", "valentine"]}, {"title": "fas fa-heart-broken", "searchTerms": ["breakup", "crushed", "dislike", "dumped", "grief", "love", "lovesick", "relationship", "sad"]}, {"title": "fas fa-heartbeat", "searchTerms": ["ekg", "electrocardiogram", "health", "lifeline", "vital signs"]}, {"title": "fas fa-helicopter", "searchTerms": ["airwolf", "apache", "chopper", "flight", "fly", "travel"]}, {"title": "fas fa-highlighter", "searchTerms": ["edit", "marker", "sharpie", "update", "write"]}, {"title": "fas fa-hiking", "searchTerms": ["activity", "backpack", "fall", "fitness", "outdoors", "person", "seasonal", "walking"]}, {"title": "fas fa-hippo", "searchTerms": ["animal", "fauna", "hippopotamus", "hungry", "mammal"]}, {"title": "fab fa-hips", "searchTerms": []}, {"title": "fab fa-hire-a-helper", "searchTerms": []}, {"title": "fas fa-history", "searchTerms": ["Rewind", "clock", "reverse", "time", "time machine"]}, {"title": "fas fa-hockey-puck", "searchTerms": ["ice", "nhl", "sport"]}, {"title": "fas fa-holly-berry", "searchTerms": ["catwoman", "christmas", "decoration", "flora", "halle", "holiday", "ororo munroe", "plant", "storm", "xmas"]}, {"title": "fas fa-home", "searchTerms": ["abode", "building", "house", "main"]}, {"title": "fab fa-hooli", "searchTerms": []}, {"title": "fab fa-hornbill", "searchTerms": []}, {"title": "fas fa-horse", "searchTerms": ["equus", "fauna", "mammmal", "mare", "neigh", "pony"]}, {"title": "fas fa-horse-head", "searchTerms": ["equus", "fauna", "mammmal", "mare", "neigh", "pony"]}, {"title": "fas fa-hospital", "searchTerms": ["building", "emergency room", "medical center"]}, {"title": "far fa-hospital", "searchTerms": ["building", "emergency room", "medical center"]}, {"title": "fas fa-hospital-alt", "searchTerms": ["building", "emergency room", "medical center"]}, {"title": "fas fa-hospital-symbol", "searchTerms": ["clinic", "emergency", "map"]}, {"title": "fas fa-hot-tub", "searchTerms": ["bath", "jacuzzi", "massage", "sauna", "spa"]}, {"title": "fas fa-hotdog", "searchTerms": ["bun", "chili", "frankfurt", "frankfurter", "kosher", "polish", "sandwich", "sausage", "vienna", "weiner"]}, {"title": "fas fa-hotel", "searchTerms": ["building", "inn", "lodging", "motel", "resort", "travel"]}, {"title": "fab fa-hotjar", "searchTerms": []}, {"title": "fas fa-hourglass", "searchTerms": ["hour", "minute", "sand", "stopwatch", "time"]}, {"title": "far fa-hourglass", "searchTerms": ["hour", "minute", "sand", "stopwatch", "time"]}, {"title": "fas fa-hourglass-end", "searchTerms": ["hour", "minute", "sand", "stopwatch", "time"]}, {"title": "fas fa-hourglass-half", "searchTerms": ["hour", "minute", "sand", "stopwatch", "time"]}, {"title": "fas fa-hourglass-start", "searchTerms": ["hour", "minute", "sand", "stopwatch", "time"]}, {"title": "fas fa-house-damage", "searchTerms": ["building", "devastation", "disaster", "home", "insurance"]}, {"title": "fab fa-houzz", "searchTerms": []}, {"title": "fas fa-hryvnia", "searchTerms": ["currency", "money", "ukraine", "ukrainian"]}, {"title": "fab fa-html5", "searchTerms": []}, {"title": "fab fa-hubspot", "searchTerms": []}, {"title": "fas fa-i-cursor", "searchTerms": ["editing", "i-beam", "type", "writing"]}, {"title": "fas fa-ice-cream", "searchTerms": ["chocolate", "cone", "dessert", "frozen", "scoop", "sorbet", "vanilla", "yogurt"]}, {"title": "fas fa-icicles", "searchTerms": ["cold", "frozen", "hanging", "ice", "seasonal", "sharp"]}, {"title": "fas fa-icons", "searchTerms": ["bolt", "emoji", "heart", "image", "music", "photo", "symbols"]}, {"title": "fas fa-id-badge", "searchTerms": ["address", "contact", "identification", "license", "profile"]}, {"title": "far fa-id-badge", "searchTerms": ["address", "contact", "identification", "license", "profile"]}, {"title": "fas fa-id-card", "searchTerms": ["contact", "demographics", "document", "identification", "issued", "profile"]}, {"title": "far fa-id-card", "searchTerms": ["contact", "demographics", "document", "identification", "issued", "profile"]}, {"title": "fas fa-id-card-alt", "searchTerms": ["contact", "demographics", "document", "identification", "issued", "profile"]}, {"title": "fab fa-ideal", "searchTerms": []}, {"title": "fas fa-igloo", "searchTerms": ["dome", "dwelling", "eskimo", "home", "house", "ice", "snow"]}, {"title": "fas fa-image", "searchTerms": ["album", "landscape", "photo", "picture"]}, {"title": "far fa-image", "searchTerms": ["album", "landscape", "photo", "picture"]}, {"title": "fas fa-images", "searchTerms": ["album", "landscape", "photo", "picture"]}, {"title": "far fa-images", "searchTerms": ["album", "landscape", "photo", "picture"]}, {"title": "fab fa-imdb", "searchTerms": []}, {"title": "fas fa-inbox", "searchTerms": ["archive", "desk", "email", "mail", "message"]}, {"title": "fas fa-indent", "searchTerms": ["align", "justify", "paragraph", "tab"]}, {"title": "fas fa-industry", "searchTerms": ["building", "factory", "industrial", "manufacturing", "mill", "warehouse"]}, {"title": "fas fa-infinity", "searchTerms": ["eternity", "forever", "math"]}, {"title": "fas fa-info", "searchTerms": ["details", "help", "information", "more", "support"]}, {"title": "fas fa-info-circle", "searchTerms": ["details", "help", "information", "more", "support"]}, {"title": "fab fa-instagram", "searchTerms": []}, {"title": "fab fa-instagram-square", "searchTerms": []}, {"title": "fab fa-intercom", "searchTerms": ["app", "customer", "messenger"]}, {"title": "fab fa-internet-explorer", "searchTerms": ["browser", "ie"]}, {"title": "fab fa-invision", "searchTerms": ["app", "design", "interface"]}, {"title": "fab fa-ioxhost", "searchTerms": []}, {"title": "fas fa-italic", "searchTerms": ["edit", "emphasis", "font", "format", "text", "type"]}, {"title": "fab fa-itch-io", "searchTerms": []}, {"title": "fab fa-itunes", "searchTerms": []}, {"title": "fab fa-itunes-note", "searchTerms": []}, {"title": "fab fa-java", "searchTerms": []}, {"title": "fas fa-jedi", "searchTerms": ["crest", "force", "sith", "skywalker", "star wars", "yoda"]}, {"title": "fab fa-jedi-order", "searchTerms": ["star wars"]}, {"title": "fab fa-jenkins", "searchTerms": []}, {"title": "fab fa-jira", "searchTerms": ["atlassian"]}, {"title": "fab fa-joget", "searchTerms": []}, {"title": "fas fa-joint", "searchTerms": ["blunt", "cannabis", "doobie", "drugs", "marijuana", "roach", "smoke", "smoking", "spliff"]}, {"title": "fab fa-joomla", "searchTerms": []}, {"title": "fas fa-journal-whills", "searchTerms": ["book", "force", "jedi", "sith", "star wars", "yoda"]}, {"title": "fab fa-js", "searchTerms": []}, {"title": "fab fa-js-square", "searchTerms": []}, {"title": "fab fa-jsfiddle", "searchTerms": []}, {"title": "fas fa-kaaba", "searchTerms": ["building", "cube", "islam", "muslim"]}, {"title": "fab fa-kaggle", "searchTerms": []}, {"title": "fas fa-key", "searchTerms": ["lock", "password", "private", "secret", "unlock"]}, {"title": "fab fa-keybase", "searchTerms": []}, {"title": "fas fa-keyboard", "searchTerms": ["accessory", "edit", "input", "text", "type", "write"]}, {"title": "far fa-keyboard", "searchTerms": ["accessory", "edit", "input", "text", "type", "write"]}, {"title": "fab fa-keycdn", "searchTerms": []}, {"title": "fas fa-khanda", "searchTerms": ["chakkar", "sikh", "sikhism", "sword"]}, {"title": "fab fa-kickstarter", "searchTerms": []}, {"title": "fab fa-kickstarter-k", "searchTerms": []}, {"title": "fas fa-kiss", "searchTerms": ["beso", "emoticon", "face", "love", "smooch"]}, {"title": "far fa-kiss", "searchTerms": ["beso", "emoticon", "face", "love", "smooch"]}, {"title": "fas fa-kiss-beam", "searchTerms": ["beso", "emoticon", "face", "love", "smooch"]}, {"title": "far fa-kiss-beam", "searchTerms": ["beso", "emoticon", "face", "love", "smooch"]}, {"title": "fas fa-kiss-wink-heart", "searchTerms": ["beso", "emoticon", "face", "love", "smooch"]}, {"title": "far fa-kiss-wink-heart", "searchTerms": ["beso", "emoticon", "face", "love", "smooch"]}, {"title": "fas fa-kiwi-bird", "searchTerms": ["bird", "fauna", "new zealand"]}, {"title": "fab fa-korvue", "searchTerms": []}, {"title": "fas fa-landmark", "searchTerms": ["building", "historic", "memorable", "monument", "politics"]}, {"title": "fas fa-language", "searchTerms": ["dialect", "idiom", "localize", "speech", "translate", "vernacular"]}, {"title": "fas fa-laptop", "searchTerms": ["computer", "cpu", "dell", "demo", "device", "mac", "macbook", "machine", "pc"]}, {"title": "fas fa-laptop-code", "searchTerms": ["computer", "cpu", "dell", "demo", "develop", "device", "mac", "macbook", "machine", "pc"]}, {"title": "fas fa-laptop-medical", "searchTerms": ["computer", "device", "ehr", "electronic health records", "history"]}, {"title": "fab fa-laravel", "searchTerms": []}, {"title": "fab fa-lastfm", "searchTerms": []}, {"title": "fab fa-lastfm-square", "searchTerms": []}, {"title": "fas fa-laugh", "searchTerms": ["LOL", "emoticon", "face", "laugh", "smile"]}, {"title": "far fa-laugh", "searchTerms": ["LOL", "emoticon", "face", "laugh", "smile"]}, {"title": "fas fa-laugh-beam", "searchTerms": ["LOL", "emoticon", "face", "happy", "smile"]}, {"title": "far fa-laugh-beam", "searchTerms": ["LOL", "emoticon", "face", "happy", "smile"]}, {"title": "fas fa-laugh-squint", "searchTerms": ["LOL", "emoticon", "face", "happy", "smile"]}, {"title": "far fa-laugh-squint", "searchTerms": ["LOL", "emoticon", "face", "happy", "smile"]}, {"title": "fas fa-laugh-wink", "searchTerms": ["LOL", "emoticon", "face", "happy", "smile"]}, {"title": "far fa-laugh-wink", "searchTerms": ["LOL", "emoticon", "face", "happy", "smile"]}, {"title": "fas fa-layer-group", "searchTerms": ["arrange", "develop", "layers", "map", "stack"]}, {"title": "fas fa-leaf", "searchTerms": ["eco", "flora", "nature", "plant", "vegan"]}, {"title": "fab fa-leanpub", "searchTerms": []}, {"title": "fas fa-lemon", "searchTerms": ["citrus", "lemonade", "lime", "tart"]}, {"title": "far fa-lemon", "searchTerms": ["citrus", "lemonade", "lime", "tart"]}, {"title": "fab fa-less", "searchTerms": []}, {"title": "fas fa-less-than", "searchTerms": ["arithmetic", "compare", "math"]}, {"title": "fas fa-less-than-equal", "searchTerms": ["arithmetic", "compare", "math"]}, {"title": "fas fa-level-down-alt", "searchTerms": ["arrow", "level-down"]}, {"title": "fas fa-level-up-alt", "searchTerms": ["arrow", "level-up"]}, {"title": "fas fa-life-ring", "searchTerms": ["coast guard", "help", "overboard", "save", "support"]}, {"title": "far fa-life-ring", "searchTerms": ["coast guard", "help", "overboard", "save", "support"]}, {"title": "fas fa-lightbulb", "searchTerms": ["energy", "idea", "inspiration", "light"]}, {"title": "far fa-lightbulb", "searchTerms": ["energy", "idea", "inspiration", "light"]}, {"title": "fab fa-line", "searchTerms": []}, {"title": "fas fa-link", "searchTerms": ["attach", "attachment", "chain", "connect"]}, {"title": "fab fa-linkedin", "searchTerms": ["linkedin-square"]}, {"title": "fab fa-linkedin-in", "searchTerms": ["linkedin"]}, {"title": "fab fa-linode", "searchTerms": []}, {"title": "fab fa-linux", "searchTerms": ["tux"]}, {"title": "fas fa-lira-sign", "searchTerms": ["currency", "money", "try", "turkish"]}, {"title": "fas fa-list", "searchTerms": ["checklist", "completed", "done", "finished", "ol", "todo", "ul"]}, {"title": "fas fa-list-alt", "searchTerms": ["checklist", "completed", "done", "finished", "ol", "todo", "ul"]}, {"title": "far fa-list-alt", "searchTerms": ["checklist", "completed", "done", "finished", "ol", "todo", "ul"]}, {"title": "fas fa-list-ol", "searchTerms": ["checklist", "completed", "done", "finished", "numbers", "ol", "todo", "ul"]}, {"title": "fas fa-list-ul", "searchTerms": ["checklist", "completed", "done", "finished", "ol", "todo", "ul"]}, {"title": "fas fa-location-arrow", "searchTerms": ["address", "compass", "coordinate", "direction", "gps", "map", "navigation", "place"]}, {"title": "fas fa-lock", "searchTerms": ["admin", "lock", "open", "password", "private", "protect", "security"]}, {"title": "fas fa-lock-open", "searchTerms": ["admin", "lock", "open", "password", "private", "protect", "security"]}, {"title": "fas fa-long-arrow-alt-down", "searchTerms": ["download", "long-arrow-down"]}, {"title": "fas fa-long-arrow-alt-left", "searchTerms": ["back", "long-arrow-left", "previous"]}, {"title": "fas fa-long-arrow-alt-right", "searchTerms": ["forward", "long-arrow-right", "next"]}, {"title": "fas fa-long-arrow-alt-up", "searchTerms": ["long-arrow-up", "upload"]}, {"title": "fas fa-low-vision", "searchTerms": ["blind", "eye", "sight"]}, {"title": "fas fa-luggage-cart", "searchTerms": ["bag", "baggage", "suitcase", "travel"]}, {"title": "fab fa-lyft", "searchTerms": []}, {"title": "fab fa-magento", "searchTerms": []}, {"title": "fas fa-magic", "searchTerms": ["autocomplete", "automatic", "mage", "magic", "spell", "wand", "witch", "wizard"]}, {"title": "fas fa-magnet", "searchTerms": ["Attract", "lodestone", "tool"]}, {"title": "fas fa-mail-bulk", "searchTerms": ["archive", "envelope", "letter", "post office", "postal", "postcard", "send", "stamp", "usps"]}, {"title": "fab fa-mailchimp", "searchTerms": []}, {"title": "fas fa-male", "searchTerms": ["human", "man", "person", "profile", "user"]}, {"title": "fab fa-mandalorian", "searchTerms": []}, {"title": "fas fa-map", "searchTerms": ["address", "coordinates", "destination", "gps", "localize", "location", "map", "navigation", "paper", "pin", "place", "point of interest", "position", "route", "travel"]}, {"title": "far fa-map", "searchTerms": ["address", "coordinates", "destination", "gps", "localize", "location", "map", "navigation", "paper", "pin", "place", "point of interest", "position", "route", "travel"]}, {"title": "fas fa-map-marked", "searchTerms": ["address", "coordinates", "destination", "gps", "localize", "location", "map", "navigation", "paper", "pin", "place", "point of interest", "position", "route", "travel"]}, {"title": "fas fa-map-marked-alt", "searchTerms": ["address", "coordinates", "destination", "gps", "localize", "location", "map", "navigation", "paper", "pin", "place", "point of interest", "position", "route", "travel"]}, {"title": "fas fa-map-marker", "searchTerms": ["address", "coordinates", "destination", "gps", "localize", "location", "map", "navigation", "paper", "pin", "place", "point of interest", "position", "route", "travel"]}, {"title": "fas fa-map-marker-alt", "searchTerms": ["address", "coordinates", "destination", "gps", "localize", "location", "map", "navigation", "paper", "pin", "place", "point of interest", "position", "route", "travel"]}, {"title": "fas fa-map-pin", "searchTerms": ["address", "agree", "coordinates", "destination", "gps", "localize", "location", "map", "marker", "navigation", "pin", "place", "position", "travel"]}, {"title": "fas fa-map-signs", "searchTerms": ["directions", "directory", "map", "signage", "wayfinding"]}, {"title": "fab fa-markdown", "searchTerms": []}, {"title": "fas fa-marker", "searchTerms": ["design", "edit", "sharpie", "update", "write"]}, {"title": "fas fa-mars", "searchTerms": ["male"]}, {"title": "fas fa-mars-double", "searchTerms": []}, {"title": "fas fa-mars-stroke", "searchTerms": []}, {"title": "fas fa-mars-stroke-h", "searchTerms": []}, {"title": "fas fa-mars-stroke-v", "searchTerms": []}, {"title": "fas fa-mask", "searchTerms": ["carnivale", "costume", "disguise", "halloween", "secret", "super hero"]}, {"title": "fab fa-mastodon", "searchTerms": []}, {"title": "fab fa-maxcdn", "searchTerms": []}, {"title": "fab fa-mdb", "searchTerms": []}, {"title": "fas fa-medal", "searchTerms": ["award", "ribbon", "star", "trophy"]}, {"title": "fab fa-medapps", "searchTerms": []}, {"title": "fab fa-medium", "searchTerms": []}, {"title": "fab fa-medium-m", "searchTerms": []}, {"title": "fas fa-medkit", "searchTerms": ["first aid", "firstaid", "health", "help", "support"]}, {"title": "fab fa-medrt", "searchTerms": []}, {"title": "fab fa-meetup", "searchTerms": []}, {"title": "fab fa-megaport", "searchTerms": []}, {"title": "fas fa-meh", "searchTerms": ["emoticon", "face", "neutral", "rating"]}, {"title": "far fa-meh", "searchTerms": ["emoticon", "face", "neutral", "rating"]}, {"title": "fas fa-meh-blank", "searchTerms": ["emoticon", "face", "neutral", "rating"]}, {"title": "far fa-meh-blank", "searchTerms": ["emoticon", "face", "neutral", "rating"]}, {"title": "fas fa-meh-rolling-eyes", "searchTerms": ["emoticon", "face", "neutral", "rating"]}, {"title": "far fa-meh-rolling-eyes", "searchTerms": ["emoticon", "face", "neutral", "rating"]}, {"title": "fas fa-memory", "searchTerms": ["DIMM", "RAM", "hardware", "storage", "technology"]}, {"title": "fab fa-mendeley", "searchTerms": []}, {"title": "fas fa-menorah", "searchTerms": ["candle", "hanukkah", "jewish", "judaism", "light"]}, {"title": "fas fa-mercury", "searchTerms": ["transgender"]}, {"title": "fas fa-meteor", "searchTerms": ["armageddon", "asteroid", "comet", "shooting star", "space"]}, {"title": "fab fa-microblog", "searchTerms": []}, {"title": "fas fa-microchip", "searchTerms": ["cpu", "hardware", "processor", "technology"]}, {"title": "fas fa-microphone", "searchTerms": ["audio", "podcast", "record", "sing", "sound", "voice"]}, {"title": "fas fa-microphone-alt", "searchTerms": ["audio", "podcast", "record", "sing", "sound", "voice"]}, {"title": "fas fa-microphone-alt-slash", "searchTerms": ["audio", "disable", "mute", "podcast", "record", "sing", "sound", "voice"]}, {"title": "fas fa-microphone-slash", "searchTerms": ["audio", "disable", "mute", "podcast", "record", "sing", "sound", "voice"]}, {"title": "fas fa-microscope", "searchTerms": ["electron", "lens", "optics", "science", "shrink"]}, {"title": "fab fa-microsoft", "searchTerms": []}, {"title": "fas fa-minus", "searchTerms": ["collapse", "delete", "hide", "minify", "negative", "remove", "trash"]}, {"title": "fas fa-minus-circle", "searchTerms": ["delete", "hide", "negative", "remove", "shape", "trash"]}, {"title": "fas fa-minus-square", "searchTerms": ["collapse", "delete", "hide", "minify", "negative", "remove", "shape", "trash"]}, {"title": "far fa-minus-square", "searchTerms": ["collapse", "delete", "hide", "minify", "negative", "remove", "shape", "trash"]}, {"title": "fas fa-mitten", "searchTerms": ["clothing", "cold", "glove", "hands", "knitted", "seasonal", "warmth"]}, {"title": "fab fa-mix", "searchTerms": []}, {"title": "fab fa-mixcloud", "searchTerms": []}, {"title": "fab fa-mixer", "searchTerms": []}, {"title": "fab fa-mizuni", "searchTerms": []}, {"title": "fas fa-mobile", "searchTerms": ["apple", "call", "cell phone", "cellphone", "device", "iphone", "number", "screen", "telephone"]}, {"title": "fas fa-mobile-alt", "searchTerms": ["apple", "call", "cell phone", "cellphone", "device", "iphone", "number", "screen", "telephone"]}, {"title": "fab fa-modx", "searchTerms": []}, {"title": "fab fa-monero", "searchTerms": []}, {"title": "fas fa-money-bill", "searchTerms": ["buy", "cash", "checkout", "money", "payment", "price", "purchase"]}, {"title": "fas fa-money-bill-alt", "searchTerms": ["buy", "cash", "checkout", "money", "payment", "price", "purchase"]}, {"title": "far fa-money-bill-alt", "searchTerms": ["buy", "cash", "checkout", "money", "payment", "price", "purchase"]}, {"title": "fas fa-money-bill-wave", "searchTerms": ["buy", "cash", "checkout", "money", "payment", "price", "purchase"]}, {"title": "fas fa-money-bill-wave-alt", "searchTerms": ["buy", "cash", "checkout", "money", "payment", "price", "purchase"]}, {"title": "fas fa-money-check", "searchTerms": ["bank check", "buy", "checkout", "cheque", "money", "payment", "price", "purchase"]}, {"title": "fas fa-money-check-alt", "searchTerms": ["bank check", "buy", "checkout", "cheque", "money", "payment", "price", "purchase"]}, {"title": "fas fa-monument", "searchTerms": ["building", "historic", "landmark", "memorable"]}, {"title": "fas fa-moon", "searchTerms": ["contrast", "crescent", "dark", "lunar", "night"]}, {"title": "far fa-moon", "searchTerms": ["contrast", "crescent", "dark", "lunar", "night"]}, {"title": "fas fa-mortar-pestle", "searchTerms": ["crush", "culinary", "grind", "medical", "mix", "pharmacy", "prescription", "spices"]}, {"title": "fas fa-mosque", "searchTerms": ["building", "islam", "landmark", "muslim"]}, {"title": "fas fa-motorcycle", "searchTerms": ["bike", "machine", "transportation", "vehicle"]}, {"title": "fas fa-mountain", "searchTerms": ["glacier", "hiking", "hill", "landscape", "travel", "view"]}, {"title": "fas fa-mouse", "searchTerms": ["click", "computer", "cursor", "input", "peripheral"]}, {"title": "fas fa-mouse-pointer", "searchTerms": ["arrow", "cursor", "select"]}, {"title": "fas fa-mug-hot", "searchTerms": ["caliente", "cocoa", "coffee", "cup", "drink", "holiday", "hot chocolate", "steam", "tea", "warmth"]}, {"title": "fas fa-music", "searchTerms": ["lyrics", "melody", "note", "sing", "sound"]}, {"title": "fab fa-napster", "searchTerms": []}, {"title": "fab fa-neos", "searchTerms": []}, {"title": "fas fa-network-wired", "searchTerms": ["computer", "connect", "ethernet", "internet", "intranet"]}, {"title": "fas fa-neuter", "searchTerms": []}, {"title": "fas fa-newspaper", "searchTerms": ["article", "editorial", "headline", "journal", "journalism", "news", "press"]}, {"title": "far fa-newspaper", "searchTerms": ["article", "editorial", "headline", "journal", "journalism", "news", "press"]}, {"title": "fab fa-nimblr", "searchTerms": []}, {"title": "fab fa-node", "searchTerms": []}, {"title": "fab fa-node-js", "searchTerms": []}, {"title": "fas fa-not-equal", "searchTerms": ["arithmetic", "compare", "math"]}, {"title": "fas fa-notes-medical", "searchTerms": ["clipboard", "doctor", "ehr", "health", "history", "records"]}, {"title": "fab fa-npm", "searchTerms": []}, {"title": "fab fa-ns8", "searchTerms": []}, {"title": "fab fa-nutritionix", "searchTerms": []}, {"title": "fas fa-object-group", "searchTerms": ["combine", "copy", "design", "merge", "select"]}, {"title": "far fa-object-group", "searchTerms": ["combine", "copy", "design", "merge", "select"]}, {"title": "fas fa-object-ungroup", "searchTerms": ["copy", "design", "merge", "select", "separate"]}, {"title": "far fa-object-ungroup", "searchTerms": ["copy", "design", "merge", "select", "separate"]}, {"title": "fab fa-odnoklassniki", "searchTerms": []}, {"title": "fab fa-odnoklassniki-square", "searchTerms": []}, {"title": "fas fa-oil-can", "searchTerms": ["auto", "crude", "gasoline", "grease", "lubricate", "petroleum"]}, {"title": "fab fa-old-republic", "searchTerms": ["politics", "star wars"]}, {"title": "fas fa-om", "searchTerms": ["buddhism", "hinduism", "jainism", "mantra"]}, {"title": "fab fa-opencart", "searchTerms": []}, {"title": "fab fa-openid", "searchTerms": []}, {"title": "fab fa-opera", "searchTerms": []}, {"title": "fab fa-optin-monster", "searchTerms": []}, {"title": "fab fa-orcid", "searchTerms": []}, {"title": "fab fa-osi", "searchTerms": []}, {"title": "fas fa-otter", "searchTerms": ["animal", "badger", "fauna", "fur", "mammal", "marten"]}, {"title": "fas fa-outdent", "searchTerms": ["align", "justify", "paragraph", "tab"]}, {"title": "fab fa-page4", "searchTerms": []}, {"title": "fab fa-pagelines", "searchTerms": ["eco", "flora", "leaf", "leaves", "nature", "plant", "tree"]}, {"title": "fas fa-pager", "searchTerms": ["beeper", "cellphone", "communication"]}, {"title": "fas fa-paint-brush", "searchTerms": ["acrylic", "art", "brush", "color", "fill", "paint", "pigment", "watercolor"]}, {"title": "fas fa-paint-roller", "searchTerms": ["acrylic", "art", "brush", "color", "fill", "paint", "pigment", "watercolor"]}, {"title": "fas fa-palette", "searchTerms": ["acrylic", "art", "brush", "color", "fill", "paint", "pigment", "watercolor"]}, {"title": "fab fa-palfed", "searchTerms": []}, {"title": "fas fa-pallet", "searchTerms": ["archive", "box", "inventory", "shipping", "warehouse"]}, {"title": "fas fa-paper-plane", "searchTerms": ["air", "float", "fold", "mail", "paper", "send"]}, {"title": "far fa-paper-plane", "searchTerms": ["air", "float", "fold", "mail", "paper", "send"]}, {"title": "fas fa-paperclip", "searchTerms": ["attach", "attachment", "connect", "link"]}, {"title": "fas fa-parachute-box", "searchTerms": ["aid", "assistance", "rescue", "supplies"]}, {"title": "fas fa-paragraph", "searchTerms": ["edit", "format", "text", "writing"]}, {"title": "fas fa-parking", "searchTerms": ["auto", "car", "garage", "meter"]}, {"title": "fas fa-passport", "searchTerms": ["document", "id", "identification", "issued", "travel"]}, {"title": "fas fa-pastafarianism", "searchTerms": ["agnosticism", "atheism", "flying spaghetti monster", "fsm"]}, {"title": "fas fa-paste", "searchTerms": ["clipboard", "copy", "document", "paper"]}, {"title": "fab fa-patreon", "searchTerms": []}, {"title": "fas fa-pause", "searchTerms": ["hold", "wait"]}, {"title": "fas fa-pause-circle", "searchTerms": ["hold", "wait"]}, {"title": "far fa-pause-circle", "searchTerms": ["hold", "wait"]}, {"title": "fas fa-paw", "searchTerms": ["animal", "cat", "dog", "pet", "print"]}, {"title": "fab fa-paypal", "searchTerms": []}, {"title": "fas fa-peace", "searchTerms": ["serenity", "tranquility", "truce", "war"]}, {"title": "fas fa-pen", "searchTerms": ["design", "edit", "update", "write"]}, {"title": "fas fa-pen-alt", "searchTerms": ["design", "edit", "update", "write"]}, {"title": "fas fa-pen-fancy", "searchTerms": ["design", "edit", "fountain pen", "update", "write"]}, {"title": "fas fa-pen-nib", "searchTerms": ["design", "edit", "fountain pen", "update", "write"]}, {"title": "fas fa-pen-square", "searchTerms": ["edit", "pencil-square", "update", "write"]}, {"title": "fas fa-pencil-alt", "searchTerms": ["design", "edit", "pencil", "update", "write"]}, {"title": "fas fa-pencil-ruler", "searchTerms": ["design", "draft", "draw", "pencil"]}, {"title": "fab fa-penny-arcade", "searchTerms": ["Dungeons & Dragons", "d&d", "dnd", "fantasy", "game", "gaming", "pax", "tabletop"]}, {"title": "fas fa-people-carry", "searchTerms": ["box", "carry", "fragile", "help", "movers", "package"]}, {"title": "fas fa-pepper-hot", "searchTerms": ["buffalo wings", "capsicum", "chili", "chilli", "habanero", "jalapeno", "mexican", "spicy", "tabasco", "vegetable"]}, {"title": "fas fa-percent", "searchTerms": ["discount", "fraction", "proportion", "rate", "ratio"]}, {"title": "fas fa-percentage", "searchTerms": ["discount", "fraction", "proportion", "rate", "ratio"]}, {"title": "fab fa-periscope", "searchTerms": []}, {"title": "fas fa-person-booth", "searchTerms": ["changing", "changing room", "election", "human", "person", "vote", "voting"]}, {"title": "fab fa-phabricator", "searchTerms": []}, {"title": "fab fa-phoenix-framework", "searchTerms": []}, {"title": "fab fa-phoenix-squadron", "searchTerms": []}, {"title": "fas fa-phone", "searchTerms": ["call", "earphone", "number", "support", "telephone", "voice"]}, {"title": "fas fa-phone-alt", "searchTerms": ["call", "earphone", "number", "support", "telephone", "voice"]}, {"title": "fas fa-phone-slash", "searchTerms": ["call", "cancel", "earphone", "mute", "number", "support", "telephone", "voice"]}, {"title": "fas fa-phone-square", "searchTerms": ["call", "earphone", "number", "support", "telephone", "voice"]}, {"title": "fas fa-phone-square-alt", "searchTerms": ["call", "earphone", "number", "support", "telephone", "voice"]}, {"title": "fas fa-phone-volume", "searchTerms": ["call", "earphone", "number", "sound", "support", "telephone", "voice", "volume-control-phone"]}, {"title": "fas fa-photo-video", "searchTerms": ["av", "film", "image", "library", "media"]}, {"title": "fab fa-php", "searchTerms": []}, {"title": "fab fa-pied-piper", "searchTerms": []}, {"title": "fab fa-pied-piper-alt", "searchTerms": []}, {"title": "fab fa-pied-piper-hat", "searchTerms": ["clothing"]}, {"title": "fab fa-pied-piper-pp", "searchTerms": []}, {"title": "fab fa-pied-piper-square", "searchTerms": []}, {"title": "fas fa-piggy-bank", "searchTerms": ["bank", "save", "savings"]}, {"title": "fas fa-pills", "searchTerms": ["drugs", "medicine", "prescription", "tablets"]}, {"title": "fab fa-pinterest", "searchTerms": []}, {"title": "fab fa-pinterest-p", "searchTerms": []}, {"title": "fab fa-pinterest-square", "searchTerms": []}, {"title": "fas fa-pizza-slice", "searchTerms": ["cheese", "chicago", "italian", "mozzarella", "new york", "pepperoni", "pie", "slice", "teenage mutant ninja turtles", "tomato"]}, {"title": "fas fa-place-of-worship", "searchTerms": ["building", "church", "holy", "mosque", "synagogue"]}, {"title": "fas fa-plane", "searchTerms": ["airplane", "destination", "fly", "location", "mode", "travel", "trip"]}, {"title": "fas fa-plane-arrival", "searchTerms": ["airplane", "arriving", "destination", "fly", "land", "landing", "location", "mode", "travel", "trip"]}, {"title": "fas fa-plane-departure", "searchTerms": ["airplane", "departing", "destination", "fly", "location", "mode", "take off", "taking off", "travel", "trip"]}, {"title": "fas fa-play", "searchTerms": ["audio", "music", "playing", "sound", "start", "video"]}, {"title": "fas fa-play-circle", "searchTerms": ["audio", "music", "playing", "sound", "start", "video"]}, {"title": "far fa-play-circle", "searchTerms": ["audio", "music", "playing", "sound", "start", "video"]}, {"title": "fab fa-playstation", "searchTerms": []}, {"title": "fas fa-plug", "searchTerms": ["connect", "electric", "online", "power"]}, {"title": "fas fa-plus", "searchTerms": ["add", "create", "expand", "new", "positive", "shape"]}, {"title": "fas fa-plus-circle", "searchTerms": ["add", "create", "expand", "new", "positive", "shape"]}, {"title": "fas fa-plus-square", "searchTerms": ["add", "create", "expand", "new", "positive", "shape"]}, {"title": "far fa-plus-square", "searchTerms": ["add", "create", "expand", "new", "positive", "shape"]}, {"title": "fas fa-podcast", "searchTerms": ["audio", "broadcast", "music", "sound"]}, {"title": "fas fa-poll", "searchTerms": ["results", "survey", "trend", "vote", "voting"]}, {"title": "fas fa-poll-h", "searchTerms": ["results", "survey", "trend", "vote", "voting"]}, {"title": "fas fa-poo", "searchTerms": ["crap", "poop", "shit", "smile", "turd"]}, {"title": "fas fa-poo-storm", "searchTerms": ["bolt", "cloud", "euphemism", "lightning", "mess", "poop", "shit", "turd"]}, {"title": "fas fa-poop", "searchTerms": ["crap", "poop", "shit", "smile", "turd"]}, {"title": "fas fa-portrait", "searchTerms": ["id", "image", "photo", "picture", "selfie"]}, {"title": "fas fa-pound-sign", "searchTerms": ["currency", "gbp", "money"]}, {"title": "fas fa-power-off", "searchTerms": ["cancel", "computer", "on", "reboot", "restart"]}, {"title": "fas fa-pray", "searchTerms": ["kneel", "preach", "religion", "worship"]}, {"title": "fas fa-praying-hands", "searchTerms": ["kneel", "preach", "religion", "worship"]}, {"title": "fas fa-prescription", "searchTerms": ["drugs", "medical", "medicine", "pharmacy", "rx"]}, {"title": "fas fa-prescription-bottle", "searchTerms": ["drugs", "medical", "medicine", "pharmacy", "rx"]}, {"title": "fas fa-prescription-bottle-alt", "searchTerms": ["drugs", "medical", "medicine", "pharmacy", "rx"]}, {"title": "fas fa-print", "searchTerms": ["business", "copy", "document", "office", "paper"]}, {"title": "fas fa-procedures", "searchTerms": ["EKG", "bed", "electrocardiogram", "health", "hospital", "life", "patient", "vital"]}, {"title": "fab fa-product-hunt", "searchTerms": []}, {"title": "fas fa-project-diagram", "searchTerms": ["chart", "graph", "network", "pert"]}, {"title": "fab fa-pushed", "searchTerms": []}, {"title": "fas fa-puzzle-piece", "searchTerms": ["add-on", "addon", "game", "section"]}, {"title": "fab fa-python", "searchTerms": []}, {"title": "fab fa-qq", "searchTerms": []}, {"title": "fas fa-qrcode", "searchTerms": ["barcode", "info", "information", "scan"]}, {"title": "fas fa-question", "searchTerms": ["help", "information", "support", "unknown"]}, {"title": "fas fa-question-circle", "searchTerms": ["help", "information", "support", "unknown"]}, {"title": "far fa-question-circle", "searchTerms": ["help", "information", "support", "unknown"]}, {"title": "fas fa-quidditch", "searchTerms": ["ball", "bludger", "broom", "golden snitch", "harry potter", "hogwarts", "quaffle", "sport", "wizard"]}, {"title": "fab fa-quinscape", "searchTerms": []}, {"title": "fab fa-quora", "searchTerms": []}, {"title": "fas fa-quote-left", "searchTerms": ["mention", "note", "phrase", "text", "type"]}, {"title": "fas fa-quote-right", "searchTerms": ["mention", "note", "phrase", "text", "type"]}, {"title": "fas fa-quran", "searchTerms": ["book", "islam", "muslim", "religion"]}, {"title": "fab fa-r-project", "searchTerms": []}, {"title": "fas fa-radiation", "searchTerms": ["danger", "dangerous", "deadly", "hazard", "nuclear", "radioactive", "warning"]}, {"title": "fas fa-radiation-alt", "searchTerms": ["danger", "dangerous", "deadly", "hazard", "nuclear", "radioactive", "warning"]}, {"title": "fas fa-rainbow", "searchTerms": ["gold", "leprechaun", "prism", "rain", "sky"]}, {"title": "fas fa-random", "searchTerms": ["arrows", "shuffle", "sort", "swap", "switch", "transfer"]}, {"title": "fab fa-raspberry-pi", "searchTerms": []}, {"title": "fab fa-ravelry", "searchTerms": []}, {"title": "fab fa-react", "searchTerms": []}, {"title": "fab fa-reacteurope", "searchTerms": []}, {"title": "fab fa-readme", "searchTerms": []}, {"title": "fab fa-rebel", "searchTerms": []}, {"title": "fas fa-receipt", "searchTerms": ["check", "invoice", "money", "pay", "table"]}, {"title": "fas fa-record-vinyl", "searchTerms": ["LP", "album", "analog", "music", "phonograph", "sound"]}, {"title": "fas fa-recycle", "searchTerms": ["Waste", "compost", "garbage", "reuse", "trash"]}, {"title": "fab fa-red-river", "searchTerms": []}, {"title": "fab fa-reddit", "searchTerms": []}, {"title": "fab fa-reddit-alien", "searchTerms": []}, {"title": "fab fa-reddit-square", "searchTerms": []}, {"title": "fab fa-redhat", "searchTerms": ["linux", "operating system", "os"]}, {"title": "fas fa-redo", "searchTerms": ["forward", "refresh", "reload", "repeat"]}, {"title": "fas fa-redo-alt", "searchTerms": ["forward", "refresh", "reload", "repeat"]}, {"title": "fas fa-registered", "searchTerms": ["copyright", "mark", "trademark"]}, {"title": "far fa-registered", "searchTerms": ["copyright", "mark", "trademark"]}, {"title": "fas fa-remove-format", "searchTerms": ["cancel", "font", "format", "remove", "style", "text"]}, {"title": "fab fa-renren", "searchTerms": []}, {"title": "fas fa-reply", "searchTerms": ["mail", "message", "respond"]}, {"title": "fas fa-reply-all", "searchTerms": ["mail", "message", "respond"]}, {"title": "fab fa-replyd", "searchTerms": []}, {"title": "fas fa-republican", "searchTerms": ["american", "conservative", "election", "elephant", "politics", "republican party", "right", "right-wing", "usa"]}, {"title": "fab fa-researchgate", "searchTerms": []}, {"title": "fab fa-resolving", "searchTerms": []}, {"title": "fas fa-restroom", "searchTerms": ["bathroom", "john", "loo", "potty", "washroom", "waste", "wc"]}, {"title": "fas fa-retweet", "searchTerms": ["refresh", "reload", "share", "swap"]}, {"title": "fab fa-rev", "searchTerms": []}, {"title": "fas fa-ribbon", "searchTerms": ["badge", "cause", "lapel", "pin"]}, {"title": "fas fa-ring", "searchTerms": ["Dungeons & Dragons", "Gollum", "band", "binding", "d&d", "dnd", "engagement", "fantasy", "gold", "jewelry", "marriage", "precious"]}, {"title": "fas fa-road", "searchTerms": ["highway", "map", "pavement", "route", "street", "travel"]}, {"title": "fas fa-robot", "searchTerms": ["android", "automate", "computer", "cyborg"]}, {"title": "fas fa-rocket", "searchTerms": ["aircraft", "app", "jet", "launch", "nasa", "space"]}, {"title": "fab fa-rocketchat", "searchTerms": []}, {"title": "fab fa-rockrms", "searchTerms": []}, {"title": "fas fa-route", "searchTerms": ["directions", "navigation", "travel"]}, {"title": "fas fa-rss", "searchTerms": ["blog", "feed", "journal", "news", "writing"]}, {"title": "fas fa-rss-square", "searchTerms": ["blog", "feed", "journal", "news", "writing"]}, {"title": "fas fa-ruble-sign", "searchTerms": ["currency", "money", "rub"]}, {"title": "fas fa-ruler", "searchTerms": ["design", "draft", "length", "measure", "planning"]}, {"title": "fas fa-ruler-combined", "searchTerms": ["design", "draft", "length", "measure", "planning"]}, {"title": "fas fa-ruler-horizontal", "searchTerms": ["design", "draft", "length", "measure", "planning"]}, {"title": "fas fa-ruler-vertical", "searchTerms": ["design", "draft", "length", "measure", "planning"]}, {"title": "fas fa-running", "searchTerms": ["exercise", "health", "jog", "person", "run", "sport", "sprint"]}, {"title": "fas fa-rupee-sign", "searchTerms": ["currency", "indian", "inr", "money"]}, {"title": "fas fa-sad-cry", "searchTerms": ["emoticon", "face", "tear", "tears"]}, {"title": "far fa-sad-cry", "searchTerms": ["emoticon", "face", "tear", "tears"]}, {"title": "fas fa-sad-tear", "searchTerms": ["emoticon", "face", "tear", "tears"]}, {"title": "far fa-sad-tear", "searchTerms": ["emoticon", "face", "tear", "tears"]}, {"title": "fab fa-safari", "searchTerms": ["browser"]}, {"title": "fab fa-salesforce", "searchTerms": []}, {"title": "fab fa-sass", "searchTerms": []}, {"title": "fas fa-satellite", "searchTerms": ["communications", "hardware", "orbit", "space"]}, {"title": "fas fa-satellite-dish", "searchTerms": ["SETI", "communications", "hardware", "receiver", "saucer", "signal", "space"]}, {"title": "fas fa-save", "searchTerms": ["disk", "download", "floppy", "floppy-o"]}, {"title": "far fa-save", "searchTerms": ["disk", "download", "floppy", "floppy-o"]}, {"title": "fab fa-schlix", "searchTerms": []}, {"title": "fas fa-school", "searchTerms": ["building", "education", "learn", "student", "teacher"]}, {"title": "fas fa-screwdriver", "searchTerms": ["admin", "fix", "mechanic", "repair", "settings", "tool"]}, {"title": "fab fa-scribd", "searchTerms": []}, {"title": "fas fa-scroll", "searchTerms": ["Dungeons & Dragons", "announcement", "d&d", "dnd", "fantasy", "paper", "script"]}, {"title": "fas fa-sd-card", "searchTerms": ["image", "memory", "photo", "save"]}, {"title": "fas fa-search", "searchTerms": ["bigger", "enlarge", "find", "magnify", "preview", "zoom"]}, {"title": "fas fa-search-dollar", "searchTerms": ["bigger", "enlarge", "find", "magnify", "money", "preview", "zoom"]}, {"title": "fas fa-search-location", "searchTerms": ["bigger", "enlarge", "find", "magnify", "preview", "zoom"]}, {"title": "fas fa-search-minus", "searchTerms": ["minify", "negative", "smaller", "zoom", "zoom out"]}, {"title": "fas fa-search-plus", "searchTerms": ["bigger", "enlarge", "magnify", "positive", "zoom", "zoom in"]}, {"title": "fab fa-searchengin", "searchTerms": []}, {"title": "fas fa-seedling", "searchTerms": ["flora", "grow", "plant", "vegan"]}, {"title": "fab fa-sellcast", "searchTerms": ["eercast"]}, {"title": "fab fa-sellsy", "searchTerms": []}, {"title": "fas fa-server", "searchTerms": ["computer", "cpu", "database", "hardware", "network"]}, {"title": "fab fa-servicestack", "searchTerms": []}, {"title": "fas fa-shapes", "searchTerms": ["blocks", "build", "circle", "square", "triangle"]}, {"title": "fas fa-share", "searchTerms": ["forward", "save", "send", "social"]}, {"title": "fas fa-share-alt", "searchTerms": ["forward", "save", "send", "social"]}, {"title": "fas fa-share-alt-square", "searchTerms": ["forward", "save", "send", "social"]}, {"title": "fas fa-share-square", "searchTerms": ["forward", "save", "send", "social"]}, {"title": "far fa-share-square", "searchTerms": ["forward", "save", "send", "social"]}, {"title": "fas fa-shekel-sign", "searchTerms": ["currency", "ils", "money"]}, {"title": "fas fa-shield-alt", "searchTerms": ["achievement", "award", "block", "defend", "security", "winner"]}, {"title": "fas fa-ship", "searchTerms": ["boat", "sea", "water"]}, {"title": "fas fa-shipping-fast", "searchTerms": ["express", "fedex", "mail", "overnight", "package", "ups"]}, {"title": "fab fa-shirtsinbulk", "searchTerms": []}, {"title": "fas fa-shoe-prints", "searchTerms": ["feet", "footprints", "steps", "walk"]}, {"title": "fab fa-shopify", "searchTerms": []}, {"title": "fas fa-shopping-bag", "searchTerms": ["buy", "checkout", "grocery", "payment", "purchase"]}, {"title": "fas fa-shopping-basket", "searchTerms": ["buy", "checkout", "grocery", "payment", "purchase"]}, {"title": "fas fa-shopping-cart", "searchTerms": ["buy", "checkout", "grocery", "payment", "purchase"]}, {"title": "fab fa-shopware", "searchTerms": []}, {"title": "fas fa-shower", "searchTerms": ["bath", "clean", "faucet", "water"]}, {"title": "fas fa-shuttle-van", "searchTerms": ["airport", "machine", "public-transportation", "transportation", "travel", "vehicle"]}, {"title": "fas fa-sign", "searchTerms": ["directions", "real estate", "signage", "wayfinding"]}, {"title": "fas fa-sign-in-alt", "searchTerms": ["arrow", "enter", "join", "log in", "login", "sign in", "sign up", "sign-in", "signin", "signup"]}, {"title": "fas fa-sign-language", "searchTerms": ["Translate", "asl", "deaf", "hands"]}, {"title": "fas fa-sign-out-alt", "searchTerms": ["arrow", "exit", "leave", "log out", "logout", "sign-out"]}, {"title": "fas fa-signal", "searchTerms": ["bars", "graph", "online", "reception", "status"]}, {"title": "fas fa-signature", "searchTerms": ["John Hancock", "cursive", "name", "writing"]}, {"title": "fas fa-sim-card", "searchTerms": ["hard drive", "hardware", "portable", "storage", "technology", "tiny"]}, {"title": "fab fa-simplybuilt", "searchTerms": []}, {"title": "fab fa-sistrix", "searchTerms": []}, {"title": "fas fa-sitemap", "searchTerms": ["directory", "hierarchy", "ia", "information architecture", "organization"]}, {"title": "fab fa-sith", "searchTerms": []}, {"title": "fas fa-skating", "searchTerms": ["activity", "figure skating", "fitness", "ice", "person", "winter"]}, {"title": "fab fa-sketch", "searchTerms": ["app", "design", "interface"]}, {"title": "fas fa-skiing", "searchTerms": ["activity", "downhill", "fast", "fitness", "olympics", "outdoors", "person", "seasonal", "slalom"]}, {"title": "fas fa-skiing-nordic", "searchTerms": ["activity", "cross country", "fitness", "outdoors", "person", "seasonal"]}, {"title": "fas fa-skull", "searchTerms": ["bones", "skeleton", "x-ray", "yorick"]}, {"title": "fas fa-skull-crossbones", "searchTerms": ["Dungeons & Dragons", "alert", "bones", "d&d", "danger", "dead", "deadly", "death", "dnd", "fantasy", "halloween", "holiday", "jolly-roger", "pirate", "poison", "skeleton", "warning"]}, {"title": "fab fa-skyatlas", "searchTerms": []}, {"title": "fab fa-skype", "searchTerms": []}, {"title": "fab fa-slack", "searchTerms": ["anchor", "hash", "hashtag"]}, {"title": "fab fa-slack-hash", "searchTerms": ["anchor", "hash", "hashtag"]}, {"title": "fas fa-slash", "searchTerms": ["cancel", "close", "mute", "off", "stop", "x"]}, {"title": "fas fa-sleigh", "searchTerms": ["christmas", "claus", "fly", "holiday", "santa", "sled", "snow", "xmas"]}, {"title": "fas fa-sliders-h", "searchTerms": ["adjust", "settings", "sliders", "toggle"]}, {"title": "fab fa-slideshare", "searchTerms": []}, {"title": "fas fa-smile", "searchTerms": ["approve", "emoticon", "face", "happy", "rating", "satisfied"]}, {"title": "far fa-smile", "searchTerms": ["approve", "emoticon", "face", "happy", "rating", "satisfied"]}, {"title": "fas fa-smile-beam", "searchTerms": ["emoticon", "face", "happy", "positive"]}, {"title": "far fa-smile-beam", "searchTerms": ["emoticon", "face", "happy", "positive"]}, {"title": "fas fa-smile-wink", "searchTerms": ["emoticon", "face", "happy", "hint", "joke"]}, {"title": "far fa-smile-wink", "searchTerms": ["emoticon", "face", "happy", "hint", "joke"]}, {"title": "fas fa-smog", "searchTerms": ["dragon", "fog", "haze", "pollution", "smoke", "weather"]}, {"title": "fas fa-smoking", "searchTerms": ["cancer", "cigarette", "nicotine", "smoking status", "tobacco"]}, {"title": "fas fa-smoking-ban", "searchTerms": ["ban", "cancel", "no smoking", "non-smoking"]}, {"title": "fas fa-sms", "searchTerms": ["chat", "conversation", "message", "mobile", "notification", "phone", "sms", "texting"]}, {"title": "fab fa-snapchat", "searchTerms": []}, {"title": "fab fa-snapchat-ghost", "searchTerms": []}, {"title": "fab fa-snapchat-square", "searchTerms": []}, {"title": "fas fa-snowboarding", "searchTerms": ["activity", "fitness", "olympics", "outdoors", "person"]}, {"title": "fas fa-snowflake", "searchTerms": ["precipitation", "rain", "winter"]}, {"title": "far fa-snowflake", "searchTerms": ["precipitation", "rain", "winter"]}, {"title": "fas fa-snowman", "searchTerms": ["decoration", "frost", "frosty", "holiday"]}, {"title": "fas fa-snowplow", "searchTerms": ["clean up", "cold", "road", "storm", "winter"]}, {"title": "fas fa-socks", "searchTerms": ["business socks", "business time", "clothing", "feet", "flight of the conchords", "wednesday"]}, {"title": "fas fa-solar-panel", "searchTerms": ["clean", "eco-friendly", "energy", "green", "sun"]}, {"title": "fas fa-sort", "searchTerms": ["filter", "order"]}, {"title": "fas fa-sort-alpha-down", "searchTerms": ["alphabetical", "arrange", "filter", "order", "sort-alpha-asc"]}, {"title": "fas fa-sort-alpha-down-alt", "searchTerms": ["alphabetical", "arrange", "filter", "order", "sort-alpha-asc"]}, {"title": "fas fa-sort-alpha-up", "searchTerms": ["alphabetical", "arrange", "filter", "order", "sort-alpha-desc"]}, {"title": "fas fa-sort-alpha-up-alt", "searchTerms": ["alphabetical", "arrange", "filter", "order", "sort-alpha-desc"]}, {"title": "fas fa-sort-amount-down", "searchTerms": ["arrange", "filter", "number", "order", "sort-amount-asc"]}, {"title": "fas fa-sort-amount-down-alt", "searchTerms": ["arrange", "filter", "order", "sort-amount-asc"]}, {"title": "fas fa-sort-amount-up", "searchTerms": ["arrange", "filter", "order", "sort-amount-desc"]}, {"title": "fas fa-sort-amount-up-alt", "searchTerms": ["arrange", "filter", "order", "sort-amount-desc"]}, {"title": "fas fa-sort-down", "searchTerms": ["arrow", "descending", "filter", "order", "sort-desc"]}, {"title": "fas fa-sort-numeric-down", "searchTerms": ["arrange", "filter", "numbers", "order", "sort-numeric-asc"]}, {"title": "fas fa-sort-numeric-down-alt", "searchTerms": ["arrange", "filter", "numbers", "order", "sort-numeric-asc"]}, {"title": "fas fa-sort-numeric-up", "searchTerms": ["arrange", "filter", "numbers", "order", "sort-numeric-desc"]}, {"title": "fas fa-sort-numeric-up-alt", "searchTerms": ["arrange", "filter", "numbers", "order", "sort-numeric-desc"]}, {"title": "fas fa-sort-up", "searchTerms": ["arrow", "ascending", "filter", "order", "sort-asc"]}, {"title": "fab fa-soundcloud", "searchTerms": []}, {"title": "fab fa-sourcetree", "searchTerms": []}, {"title": "fas fa-spa", "searchTerms": ["flora", "massage", "mindfulness", "plant", "wellness"]}, {"title": "fas fa-space-shuttle", "searchTerms": ["astronaut", "machine", "nasa", "rocket", "space", "transportation"]}, {"title": "fab fa-speakap", "searchTerms": []}, {"title": "fab fa-speaker-deck", "searchTerms": []}, {"title": "fas fa-spell-check", "searchTerms": ["dictionary", "edit", "editor", "grammar", "text"]}, {"title": "fas fa-spider", "searchTerms": ["arachnid", "bug", "charlotte", "crawl", "eight", "halloween"]}, {"title": "fas fa-spinner", "searchTerms": ["circle", "loading", "progress"]}, {"title": "fas fa-splotch", "searchTerms": ["Ink", "blob", "blotch", "glob", "stain"]}, {"title": "fab fa-spotify", "searchTerms": []}, {"title": "fas fa-spray-can", "searchTerms": ["Paint", "aerosol", "design", "graffiti", "tag"]}, {"title": "fas fa-square", "searchTerms": ["block", "box", "shape"]}, {"title": "far fa-square", "searchTerms": ["block", "box", "shape"]}, {"title": "fas fa-square-full", "searchTerms": ["block", "box", "shape"]}, {"title": "fas fa-square-root-alt", "searchTerms": ["arithmetic", "calculus", "division", "math"]}, {"title": "fab fa-squarespace", "searchTerms": []}, {"title": "fab fa-stack-exchange", "searchTerms": []}, {"title": "fab fa-stack-overflow", "searchTerms": []}, {"title": "fab fa-stackpath", "searchTerms": []}, {"title": "fas fa-stamp", "searchTerms": ["art", "certificate", "imprint", "rubber", "seal"]}, {"title": "fas fa-star", "searchTerms": ["achievement", "award", "favorite", "important", "night", "rating", "score"]}, {"title": "far fa-star", "searchTerms": ["achievement", "award", "favorite", "important", "night", "rating", "score"]}, {"title": "fas fa-star-and-crescent", "searchTerms": ["islam", "muslim", "religion"]}, {"title": "fas fa-star-half", "searchTerms": ["achievement", "award", "rating", "score", "star-half-empty", "star-half-full"]}, {"title": "far fa-star-half", "searchTerms": ["achievement", "award", "rating", "score", "star-half-empty", "star-half-full"]}, {"title": "fas fa-star-half-alt", "searchTerms": ["achievement", "award", "rating", "score", "star-half-empty", "star-half-full"]}, {"title": "fas fa-star-of-david", "searchTerms": ["jewish", "judaism", "religion"]}, {"title": "fas fa-star-of-life", "searchTerms": ["doctor", "emt", "first aid", "health", "medical"]}, {"title": "fab fa-staylinked", "searchTerms": []}, {"title": "fab fa-steam", "searchTerms": []}, {"title": "fab fa-steam-square", "searchTerms": []}, {"title": "fab fa-steam-symbol", "searchTerms": []}, {"title": "fas fa-step-backward", "searchTerms": ["beginning", "first", "previous", "rewind", "start"]}, {"title": "fas fa-step-forward", "searchTerms": ["end", "last", "next"]}, {"title": "fas fa-stethoscope", "searchTerms": ["diagnosis", "doctor", "general practitioner", "hospital", "infirmary", "medicine", "office", "outpatient"]}, {"title": "fab fa-sticker-mule", "searchTerms": []}, {"title": "fas fa-sticky-note", "searchTerms": ["message", "note", "paper", "reminder", "sticker"]}, {"title": "far fa-sticky-note", "searchTerms": ["message", "note", "paper", "reminder", "sticker"]}, {"title": "fas fa-stop", "searchTerms": ["block", "box", "square"]}, {"title": "fas fa-stop-circle", "searchTerms": ["block", "box", "circle", "square"]}, {"title": "far fa-stop-circle", "searchTerms": ["block", "box", "circle", "square"]}, {"title": "fas fa-stopwatch", "searchTerms": ["clock", "reminder", "time"]}, {"title": "fas fa-store", "searchTerms": ["building", "buy", "purchase", "shopping"]}, {"title": "fas fa-store-alt", "searchTerms": ["building", "buy", "purchase", "shopping"]}, {"title": "fab fa-strava", "searchTerms": []}, {"title": "fas fa-stream", "searchTerms": ["flow", "list", "timeline"]}, {"title": "fas fa-street-view", "searchTerms": ["directions", "location", "map", "navigation"]}, {"title": "fas fa-strikethrough", "searchTerms": ["cancel", "edit", "font", "format", "text", "type"]}, {"title": "fab fa-stripe", "searchTerms": []}, {"title": "fab fa-stripe-s", "searchTerms": []}, {"title": "fas fa-stroopwafel", "searchTerms": ["caramel", "cookie", "dessert", "sweets", "waffle"]}, {"title": "fab fa-studiovinari", "searchTerms": []}, {"title": "fab fa-stumbleupon", "searchTerms": []}, {"title": "fab fa-stumbleupon-circle", "searchTerms": []}, {"title": "fas fa-subscript", "searchTerms": ["edit", "font", "format", "text", "type"]}, {"title": "fas fa-subway", "searchTerms": ["machine", "railway", "train", "transportation", "vehicle"]}, {"title": "fas fa-suitcase", "searchTerms": ["baggage", "luggage", "move", "suitcase", "travel", "trip"]}, {"title": "fas fa-suitcase-rolling", "searchTerms": ["baggage", "luggage", "move", "suitcase", "travel", "trip"]}, {"title": "fas fa-sun", "searchTerms": ["brighten", "contrast", "day", "lighter", "sol", "solar", "star", "weather"]}, {"title": "far fa-sun", "searchTerms": ["brighten", "contrast", "day", "lighter", "sol", "solar", "star", "weather"]}, {"title": "fab fa-superpowers", "searchTerms": []}, {"title": "fas fa-superscript", "searchTerms": ["edit", "exponential", "font", "format", "text", "type"]}, {"title": "fab fa-supple", "searchTerms": []}, {"title": "fas fa-surprise", "searchTerms": ["emoticon", "face", "shocked"]}, {"title": "far fa-surprise", "searchTerms": ["emoticon", "face", "shocked"]}, {"title": "fab fa-suse", "searchTerms": ["linux", "operating system", "os"]}, {"title": "fas fa-swatchbook", "searchTerms": ["Pantone", "color", "design", "hue", "palette"]}, {"title": "fab fa-swift", "searchTerms": []}, {"title": "fas fa-swimmer", "searchTerms": ["athlete", "head", "man", "olympics", "person", "pool", "water"]}, {"title": "fas fa-swimming-pool", "searchTerms": ["ladder", "recreation", "swim", "water"]}, {"title": "fab fa-symfony", "searchTerms": []}, {"title": "fas fa-synagogue", "searchTerms": ["building", "jewish", "judaism", "religion", "star of david", "temple"]}, {"title": "fas fa-sync", "searchTerms": ["exchange", "refresh", "reload", "rotate", "swap"]}, {"title": "fas fa-sync-alt", "searchTerms": ["exchange", "refresh", "reload", "rotate", "swap"]}, {"title": "fas fa-syringe", "searchTerms": ["doctor", "immunizations", "medical", "needle"]}, {"title": "fas fa-table", "searchTerms": ["data", "excel", "spreadsheet"]}, {"title": "fas fa-table-tennis", "searchTerms": ["ball", "paddle", "ping pong"]}, {"title": "fas fa-tablet", "searchTerms": ["apple", "device", "ipad", "kindle", "screen"]}, {"title": "fas fa-tablet-alt", "searchTerms": ["apple", "device", "ipad", "kindle", "screen"]}, {"title": "fas fa-tablets", "searchTerms": ["drugs", "medicine", "pills", "prescription"]}, {"title": "fas fa-tachometer-alt", "searchTerms": ["dashboard", "fast", "odometer", "speed", "speedometer"]}, {"title": "fas fa-tag", "searchTerms": ["discount", "label", "price", "shopping"]}, {"title": "fas fa-tags", "searchTerms": ["discount", "label", "price", "shopping"]}, {"title": "fas fa-tape", "searchTerms": ["design", "package", "sticky"]}, {"title": "fas fa-tasks", "searchTerms": ["checklist", "downloading", "downloads", "loading", "progress", "project management", "settings", "to do"]}, {"title": "fas fa-taxi", "searchTerms": ["cab", "cabbie", "car", "car service", "lyft", "machine", "transportation", "travel", "uber", "vehicle"]}, {"title": "fab fa-teamspeak", "searchTerms": []}, {"title": "fas fa-teeth", "searchTerms": ["bite", "dental", "dentist", "gums", "mouth", "smile", "tooth"]}, {"title": "fas fa-teeth-open", "searchTerms": ["dental", "dentist", "gums bite", "mouth", "smile", "tooth"]}, {"title": "fab fa-telegram", "searchTerms": []}, {"title": "fab fa-telegram-plane", "searchTerms": []}, {"title": "fas fa-temperature-high", "searchTerms": ["cook", "mercury", "summer", "thermometer", "warm"]}, {"title": "fas fa-temperature-low", "searchTerms": ["cold", "cool", "mercury", "thermometer", "winter"]}, {"title": "fab fa-tencent-weibo", "searchTerms": []}, {"title": "fas fa-tenge", "searchTerms": ["currency", "kazakhstan", "money", "price"]}, {"title": "fas fa-terminal", "searchTerms": ["code", "command", "console", "development", "prompt"]}, {"title": "fas fa-text-height", "searchTerms": ["edit", "font", "format", "text", "type"]}, {"title": "fas fa-text-width", "searchTerms": ["edit", "font", "format", "text", "type"]}, {"title": "fas fa-th", "searchTerms": ["blocks", "boxes", "grid", "squares"]}, {"title": "fas fa-th-large", "searchTerms": ["blocks", "boxes", "grid", "squares"]}, {"title": "fas fa-th-list", "searchTerms": ["checklist", "completed", "done", "finished", "ol", "todo", "ul"]}, {"title": "fab fa-the-red-yeti", "searchTerms": []}, {"title": "fas fa-theater-masks", "searchTerms": ["comedy", "perform", "theatre", "tragedy"]}, {"title": "fab fa-themeco", "searchTerms": []}, {"title": "fab fa-themeisle", "searchTerms": []}, {"title": "fas fa-thermometer", "searchTerms": ["mercury", "status", "temperature"]}, {"title": "fas fa-thermometer-empty", "searchTerms": ["cold", "mercury", "status", "temperature"]}, {"title": "fas fa-thermometer-full", "searchTerms": ["fever", "hot", "mercury", "status", "temperature"]}, {"title": "fas fa-thermometer-half", "searchTerms": ["mercury", "status", "temperature"]}, {"title": "fas fa-thermometer-quarter", "searchTerms": ["mercury", "status", "temperature"]}, {"title": "fas fa-thermometer-three-quarters", "searchTerms": ["mercury", "status", "temperature"]}, {"title": "fab fa-think-peaks", "searchTerms": []}, {"title": "fas fa-thumbs-down", "searchTerms": ["disagree", "disapprove", "dislike", "hand", "social", "thumbs-o-down"]}, {"title": "far fa-thumbs-down", "searchTerms": ["disagree", "disapprove", "dislike", "hand", "social", "thumbs-o-down"]}, {"title": "fas fa-thumbs-up", "searchTerms": ["agree", "approve", "favorite", "hand", "like", "ok", "okay", "social", "success", "thumbs-o-up", "yes", "you got it dude"]}, {"title": "far fa-thumbs-up", "searchTerms": ["agree", "approve", "favorite", "hand", "like", "ok", "okay", "social", "success", "thumbs-o-up", "yes", "you got it dude"]}, {"title": "fas fa-thumbtack", "searchTerms": ["coordinates", "location", "marker", "pin", "thumb-tack"]}, {"title": "fas fa-ticket-alt", "searchTerms": ["movie", "pass", "support", "ticket"]}, {"title": "fas fa-times", "searchTerms": ["close", "cross", "error", "exit", "incorrect", "notice", "notification", "notify", "problem", "wrong", "x"]}, {"title": "fas fa-times-circle", "searchTerms": ["close", "cross", "exit", "incorrect", "notice", "notification", "notify", "problem", "wrong", "x"]}, {"title": "far fa-times-circle", "searchTerms": ["close", "cross", "exit", "incorrect", "notice", "notification", "notify", "problem", "wrong", "x"]}, {"title": "fas fa-tint", "searchTerms": ["color", "drop", "droplet", "raindrop", "waterdrop"]}, {"title": "fas fa-tint-slash", "searchTerms": ["color", "drop", "droplet", "raindrop", "waterdrop"]}, {"title": "fas fa-tired", "searchTerms": ["angry", "emoticon", "face", "grumpy", "upset"]}, {"title": "far fa-tired", "searchTerms": ["angry", "emoticon", "face", "grumpy", "upset"]}, {"title": "fas fa-toggle-off", "searchTerms": ["switch"]}, {"title": "fas fa-toggle-on", "searchTerms": ["switch"]}, {"title": "fas fa-toilet", "searchTerms": ["bathroom", "flush", "john", "loo", "pee", "plumbing", "poop", "porcelain", "potty", "restroom", "throne", "washroom", "waste", "wc"]}, {"title": "fas fa-toilet-paper", "searchTerms": ["bathroom", "halloween", "holiday", "lavatory", "prank", "restroom", "roll"]}, {"title": "fas fa-toolbox", "searchTerms": ["admin", "container", "fix", "repair", "settings", "tools"]}, {"title": "fas fa-tools", "searchTerms": ["admin", "fix", "repair", "screwdriver", "settings", "tools", "wrench"]}, {"title": "fas fa-tooth", "searchTerms": ["bicuspid", "dental", "dentist", "molar", "mouth", "teeth"]}, {"title": "fas fa-torah", "searchTerms": ["book", "jewish", "judaism", "religion", "scroll"]}, {"title": "fas fa-torii-gate", "searchTerms": ["building", "shintoism"]}, {"title": "fas fa-tractor", "searchTerms": ["agriculture", "farm", "vehicle"]}, {"title": "fab fa-trade-federation", "searchTerms": []}, {"title": "fas fa-trademark", "searchTerms": ["copyright", "register", "symbol"]}, {"title": "fas fa-traffic-light", "searchTerms": ["direction", "road", "signal", "travel"]}, {"title": "fas fa-trailer", "searchTerms": ["carry", "haul", "moving", "travel"]}, {"title": "fas fa-train", "searchTerms": ["bullet", "commute", "locomotive", "railway", "subway"]}, {"title": "fas fa-tram", "searchTerms": ["crossing", "machine", "mountains", "seasonal", "transportation"]}, {"title": "fas fa-transgender", "searchTerms": ["intersex"]}, {"title": "fas fa-transgender-alt", "searchTerms": ["intersex"]}, {"title": "fas fa-trash", "searchTerms": ["delete", "garbage", "hide", "remove"]}, {"title": "fas fa-trash-alt", "searchTerms": ["delete", "garbage", "hide", "remove", "trash-o"]}, {"title": "far fa-trash-alt", "searchTerms": ["delete", "garbage", "hide", "remove", "trash-o"]}, {"title": "fas fa-trash-restore", "searchTerms": ["back", "control z", "oops", "undo"]}, {"title": "fas fa-trash-restore-alt", "searchTerms": ["back", "control z", "oops", "undo"]}, {"title": "fas fa-tree", "searchTerms": ["bark", "fall", "flora", "forest", "nature", "plant", "seasonal"]}, {"title": "fab fa-trello", "searchTerms": ["atlassian"]}, {"title": "fab fa-tripadvisor", "searchTerms": []}, {"title": "fas fa-trophy", "searchTerms": ["achievement", "award", "cup", "game", "winner"]}, {"title": "fas fa-truck", "searchTerms": ["cargo", "delivery", "shipping", "vehicle"]}, {"title": "fas fa-truck-loading", "searchTerms": ["box", "cargo", "delivery", "inventory", "moving", "rental", "vehicle"]}, {"title": "fas fa-truck-monster", "searchTerms": ["offroad", "vehicle", "wheel"]}, {"title": "fas fa-truck-moving", "searchTerms": ["cargo", "inventory", "rental", "vehicle"]}, {"title": "fas fa-truck-pickup", "searchTerms": ["cargo", "vehicle"]}, {"title": "fas fa-tshirt", "searchTerms": ["clothing", "fashion", "garment", "shirt"]}, {"title": "fas fa-tty", "searchTerms": ["communication", "deaf", "telephone", "teletypewriter", "text"]}, {"title": "fab fa-tumblr", "searchTerms": []}, {"title": "fab fa-tumblr-square", "searchTerms": []}, {"title": "fas fa-tv", "searchTerms": ["computer", "display", "monitor", "television"]}, {"title": "fab fa-twitch", "searchTerms": []}, {"title": "fab fa-twitter", "searchTerms": ["social network", "tweet"]}, {"title": "fab fa-twitter-square", "searchTerms": ["social network", "tweet"]}, {"title": "fab fa-typo3", "searchTerms": []}, {"title": "fab fa-uber", "searchTerms": []}, {"title": "fab fa-ubuntu", "searchTerms": ["linux", "operating system", "os"]}, {"title": "fab fa-uikit", "searchTerms": []}, {"title": "fab fa-umbraco", "searchTerms": []}, {"title": "fas fa-umbrella", "searchTerms": ["protection", "rain", "storm", "wet"]}, {"title": "fas fa-umbrella-beach", "searchTerms": ["protection", "recreation", "sand", "shade", "summer", "sun"]}, {"title": "fas fa-underline", "searchTerms": ["edit", "emphasis", "format", "text", "writing"]}, {"title": "fas fa-undo", "searchTerms": ["back", "control z", "exchange", "oops", "return", "rotate", "swap"]}, {"title": "fas fa-undo-alt", "searchTerms": ["back", "control z", "exchange", "oops", "return", "swap"]}, {"title": "fab fa-uniregistry", "searchTerms": []}, {"title": "fab fa-unity", "searchTerms": []}, {"title": "fas fa-universal-access", "searchTerms": ["accessibility", "hearing", "person", "seeing", "visual impairment"]}, {"title": "fas fa-university", "searchTerms": ["bank", "building", "college", "higher education - students", "institution"]}, {"title": "fas fa-unlink", "searchTerms": ["attachment", "chain", "chain-broken", "remove"]}, {"title": "fas fa-unlock", "searchTerms": ["admin", "lock", "password", "private", "protect"]}, {"title": "fas fa-unlock-alt", "searchTerms": ["admin", "lock", "password", "private", "protect"]}, {"title": "fab fa-untappd", "searchTerms": []}, {"title": "fas fa-upload", "searchTerms": ["hard drive", "import", "publish"]}, {"title": "fab fa-ups", "searchTerms": ["United Parcel Service", "package", "shipping"]}, {"title": "fab fa-usb", "searchTerms": []}, {"title": "fas fa-user", "searchTerms": ["account", "avatar", "head", "human", "man", "person", "profile"]}, {"title": "far fa-user", "searchTerms": ["account", "avatar", "head", "human", "man", "person", "profile"]}, {"title": "fas fa-user-alt", "searchTerms": ["account", "avatar", "head", "human", "man", "person", "profile"]}, {"title": "fas fa-user-alt-slash", "searchTerms": ["account", "avatar", "head", "human", "man", "person", "profile"]}, {"title": "fas fa-user-astronaut", "searchTerms": ["avatar", "clothing", "cosmonaut", "nasa", "space", "suit"]}, {"title": "fas fa-user-check", "searchTerms": ["accept", "check", "person", "verified"]}, {"title": "fas fa-user-circle", "searchTerms": ["account", "avatar", "head", "human", "man", "person", "profile"]}, {"title": "far fa-user-circle", "searchTerms": ["account", "avatar", "head", "human", "man", "person", "profile"]}, {"title": "fas fa-user-clock", "searchTerms": ["alert", "person", "remind", "time"]}, {"title": "fas fa-user-cog", "searchTerms": ["admin", "cog", "person", "settings"]}, {"title": "fas fa-user-edit", "searchTerms": ["edit", "pen", "pencil", "person", "update", "write"]}, {"title": "fas fa-user-friends", "searchTerms": ["group", "people", "person", "team", "users"]}, {"title": "fas fa-user-graduate", "searchTerms": ["cap", "clothing", "commencement", "gown", "graduation", "person", "student"]}, {"title": "fas fa-user-injured", "searchTerms": ["cast", "injury", "ouch", "patient", "person", "sling"]}, {"title": "fas fa-user-lock", "searchTerms": ["admin", "lock", "person", "private", "unlock"]}, {"title": "fas fa-user-md", "searchTerms": ["job", "medical", "nurse", "occupation", "physician", "profile", "surgeon"]}, {"title": "fas fa-user-minus", "searchTerms": ["delete", "negative", "remove"]}, {"title": "fas fa-user-ninja", "searchTerms": ["assassin", "avatar", "dangerous", "deadly", "sneaky"]}, {"title": "fas fa-user-nurse", "searchTerms": ["doctor", "midwife", "practitioner", "surgeon"]}, {"title": "fas fa-user-plus", "searchTerms": ["add", "avatar", "positive", "sign up", "signup", "team"]}, {"title": "fas fa-user-secret", "searchTerms": ["clothing", "coat", "hat", "incognito", "person", "privacy", "spy", "whisper"]}, {"title": "fas fa-user-shield", "searchTerms": ["admin", "person", "private", "protect", "safe"]}, {"title": "fas fa-user-slash", "searchTerms": ["ban", "delete", "remove"]}, {"title": "fas fa-user-tag", "searchTerms": ["avatar", "discount", "label", "person", "role", "special"]}, {"title": "fas fa-user-tie", "searchTerms": ["avatar", "business", "clothing", "formal", "professional", "suit"]}, {"title": "fas fa-user-times", "searchTerms": ["archive", "delete", "remove", "x"]}, {"title": "fas fa-users", "searchTerms": ["friends", "group", "people", "persons", "profiles", "team"]}, {"title": "fas fa-users-cog", "searchTerms": ["admin", "cog", "group", "person", "settings", "team"]}, {"title": "fab fa-usps", "searchTerms": ["american", "package", "shipping", "usa"]}, {"title": "fab fa-ussunnah", "searchTerms": []}, {"title": "fas fa-utensil-spoon", "searchTerms": ["cutlery", "dining", "scoop", "silverware", "spoon"]}, {"title": "fas fa-utensils", "searchTerms": ["cutlery", "dining", "dinner", "eat", "food", "fork", "knife", "restaurant"]}, {"title": "fab fa-vaadin", "searchTerms": []}, {"title": "fas fa-vector-square", "searchTerms": ["anchors", "lines", "object", "render", "shape"]}, {"title": "fas fa-venus", "searchTerms": ["female"]}, {"title": "fas fa-venus-double", "searchTerms": ["female"]}, {"title": "fas fa-venus-mars", "searchTerms": ["Gender"]}, {"title": "fab fa-viacoin", "searchTerms": []}, {"title": "fab fa-viadeo", "searchTerms": []}, {"title": "fab fa-viadeo-square", "searchTerms": []}, {"title": "fas fa-vial", "searchTerms": ["experiment", "lab", "sample", "science", "test", "test tube"]}, {"title": "fas fa-vials", "searchTerms": ["experiment", "lab", "sample", "science", "test", "test tube"]}, {"title": "fab fa-viber", "searchTerms": []}, {"title": "fas fa-video", "searchTerms": ["camera", "film", "movie", "record", "video-camera"]}, {"title": "fas fa-video-slash", "searchTerms": ["add", "create", "film", "new", "positive", "record", "video"]}, {"title": "fas fa-vihara", "searchTerms": ["buddhism", "buddhist", "building", "monastery"]}, {"title": "fab fa-vimeo", "searchTerms": []}, {"title": "fab fa-vimeo-square", "searchTerms": []}, {"title": "fab fa-vimeo-v", "searchTerms": ["vimeo"]}, {"title": "fab fa-vine", "searchTerms": []}, {"title": "fab fa-vk", "searchTerms": []}, {"title": "fab fa-vnv", "searchTerms": []}, {"title": "fas fa-voicemail", "searchTerms": ["answer", "inbox", "message", "phone"]}, {"title": "fas fa-volleyball-ball", "searchTerms": ["beach", "olympics", "sport"]}, {"title": "fas fa-volume-down", "searchTerms": ["audio", "lower", "music", "quieter", "sound", "speaker"]}, {"title": "fas fa-volume-mute", "searchTerms": ["audio", "music", "quiet", "sound", "speaker"]}, {"title": "fas fa-volume-off", "searchTerms": ["audio", "ban", "music", "mute", "quiet", "silent", "sound"]}, {"title": "fas fa-volume-up", "searchTerms": ["audio", "higher", "louder", "music", "sound", "speaker"]}, {"title": "fas fa-vote-yea", "searchTerms": ["accept", "cast", "election", "politics", "positive", "yes"]}, {"title": "fas fa-vr-cardboard", "searchTerms": ["3d", "augment", "google", "reality", "virtual"]}, {"title": "fab fa-vuejs", "searchTerms": []}, {"title": "fas fa-walking", "searchTerms": ["exercise", "health", "pedometer", "person", "steps"]}, {"title": "fas fa-wallet", "searchTerms": ["billfold", "cash", "currency", "money"]}, {"title": "fas fa-warehouse", "searchTerms": ["building", "capacity", "garage", "inventory", "storage"]}, {"title": "fas fa-water", "searchTerms": ["lake", "liquid", "ocean", "sea", "swim", "wet"]}, {"title": "fas fa-wave-square", "searchTerms": ["frequency", "pulse", "signal"]}, {"title": "fab fa-waze", "searchTerms": []}, {"title": "fab fa-weebly", "searchTerms": []}, {"title": "fab fa-weibo", "searchTerms": []}, {"title": "fas fa-weight", "searchTerms": ["health", "measurement", "scale", "weight"]}, {"title": "fas fa-weight-hanging", "searchTerms": ["anvil", "heavy", "measurement"]}, {"title": "fab fa-weixin", "searchTerms": []}, {"title": "fab fa-whatsapp", "searchTerms": []}, {"title": "fab fa-whatsapp-square", "searchTerms": []}, {"title": "fas fa-wheelchair", "searchTerms": ["accessible", "handicap", "person"]}, {"title": "fab fa-whmcs", "searchTerms": []}, {"title": "fas fa-wifi", "searchTerms": ["connection", "hotspot", "internet", "network", "wireless"]}, {"title": "fab fa-wikipedia-w", "searchTerms": []}, {"title": "fas fa-wind", "searchTerms": ["air", "blow", "breeze", "fall", "seasonal", "weather"]}, {"title": "fas fa-window-close", "searchTerms": ["browser", "cancel", "computer", "development"]}, {"title": "far fa-window-close", "searchTerms": ["browser", "cancel", "computer", "development"]}, {"title": "fas fa-window-maximize", "searchTerms": ["browser", "computer", "development", "expand"]}, {"title": "far fa-window-maximize", "searchTerms": ["browser", "computer", "development", "expand"]}, {"title": "fas fa-window-minimize", "searchTerms": ["browser", "collapse", "computer", "development"]}, {"title": "far fa-window-minimize", "searchTerms": ["browser", "collapse", "computer", "development"]}, {"title": "fas fa-window-restore", "searchTerms": ["browser", "computer", "development"]}, {"title": "far fa-window-restore", "searchTerms": ["browser", "computer", "development"]}, {"title": "fab fa-windows", "searchTerms": ["microsoft", "operating system", "os"]}, {"title": "fas fa-wine-bottle", "searchTerms": ["alcohol", "beverage", "cabernet", "drink", "glass", "grapes", "merlot", "sauvignon"]}, {"title": "fas fa-wine-glass", "searchTerms": ["alcohol", "beverage", "cabernet", "drink", "grapes", "merlot", "sauvignon"]}, {"title": "fas fa-wine-glass-alt", "searchTerms": ["alcohol", "beverage", "cabernet", "drink", "grapes", "merlot", "sauvignon"]}, {"title": "fab fa-wix", "searchTerms": []}, {"title": "fab fa-wizards-of-the-coast", "searchTerms": ["Dungeons & Dragons", "d&d", "dnd", "fantasy", "game", "gaming", "tabletop"]}, {"title": "fab fa-wolf-pack-battalion", "searchTerms": []}, {"title": "fas fa-won-sign", "searchTerms": ["currency", "krw", "money"]}, {"title": "fab fa-wordpress", "searchTerms": []}, {"title": "fab fa-wordpress-simple", "searchTerms": []}, {"title": "fab fa-wpbeginner", "searchTerms": []}, {"title": "fab fa-wpexplorer", "searchTerms": []}, {"title": "fab fa-wpforms", "searchTerms": []}, {"title": "fab fa-wpressr", "searchTerms": ["rendact"]}, {"title": "fas fa-wrench", "searchTerms": ["construction", "fix", "mechanic", "plumbing", "settings", "spanner", "tool", "update"]}, {"title": "fas fa-x-ray", "searchTerms": ["health", "medical", "radiological images", "radiology", "skeleton"]}, {"title": "fab fa-xbox", "searchTerms": []}, {"title": "fab fa-xing", "searchTerms": []}, {"title": "fab fa-xing-square", "searchTerms": []}, {"title": "fab fa-y-combinator", "searchTerms": []}, {"title": "fab fa-yahoo", "searchTerms": []}, {"title": "fab fa-yammer", "searchTerms": []}, {"title": "fab fa-yandex", "searchTerms": []}, {"title": "fab fa-yandex-international", "searchTerms": []}, {"title": "fab fa-yarn", "searchTerms": []}, {"title": "fab fa-yelp", "searchTerms": []}, {"title": "fas fa-yen-sign", "searchTerms": ["currency", "jpy", "money"]}, {"title": "fas fa-yin-yang", "searchTerms": ["daoism", "opposites", "taoism"]}, {"title": "fab fa-yoast", "searchTerms": []}, {"title": "fab fa-youtube", "searchTerms": ["film", "video", "youtube-play", "youtube-square"]}, {"title": "fab fa-youtube-square", "searchTerms": []}, {"title": "fab fa-zhihu", "searchTerms": []}]});
