/*
    自定义工具
    I am liangzhenyu
    2018-05-25
*/

/*
插件-自定义下拉选择框（美化）
兼容：IE9+
关键的：不改变开发者原有对dom节点操作习惯，在这里使用该插件只需添加一个 $.initSelect(selector) 创建后就行，然后对 原始select 原有操作依然有效；
使用方法说明：
    1.此插件基于jQuery编写，使用时需要先导入jQuery；
       * 主要是美化自带的select
    2.获取对象
        var mySelect = $.initSelect(selector, copyStyle, callback); //三个参数的排序不限制，可随意排序
            selector:css选择器，比如：.class,#id,dom
            copyStyle：是否复制 原始select 的样式，boolean类型（可选参数），默认为false
            callback:回调函数（可选参数）；回调函数中的 this 对象为 ‘mySelect’ 对象
    3.方法说明
        change(callback):当选择项改变时响应
            callback:回调函数；回调函数中的 this 对象为 ‘mySelect’ 对象
            params:回调函数中的参数,所选中的子项的属性和值
        addOption(option):添加子项（原始方式）
            option:原select下的option标签，例子："<option name='测试' value='test'>测试</option>"
        addRow(text, params):添加子项（封装方式）
            text:选项的内容，
            params:所有属性，例子：{name:'测试',value:'test'}
        getSelected(attrName):获取所选中的子项的属性值
            attrName:可选；属性名称,默认值：'value'
            return:对应属性的值
        getSelectedText():获取所选中的子项的标签内容
            return:子项内容
        getSelectedParams():获取所选中子项的封装好的属性键值对
            return:所选中的子项的属性和值，例子：{name:"xxx",abc:"xxx",def:"xxx"}
        setSelected(attrName,attrVal):设置当前选中状态
            attrName:属性名称
            attrVal:属性的值
        setSelected(attrVal):设置当前选中状态
            attrVal:属性(名称：value)的值


     注：美化下拉框和原始下拉框是相互关联的，可以不使用作者写的方法（最少要创建对象），

 */
;(function ($, window, document, undefined) {
    $.extend({
        initSelect: function (selector, copyStyle, callback) {
            var toRGB = function (color) {
                color = color.toLowerCase();
                var __EN_COLOR = function () {
                        var placeholder = "N,O,P,Q,R,S,T,U,V,W,X,Y,Z,M,L,K,J".split(",");
                        var replace = "dark,deep,light,medium,white,red,black,aqua,green,blue,yellow,orange,pink,00,ff,80,ee".split(",");
                        var strEnColor = '{"aliceW":"f0f8L","antiqueR":"faebd7","U":"MLL","Umarine":"7Lfd4","azure":"f0LL","beige":"f5f5dc","bisque":"Le4c4","T":"MMM","blanchedalmond":"Lebcd","W":"MML","Wviolet":"8a2be2","brown":"a52a2a","burlywood":"deb887","cadetW":"5f9ea0","chartreuse":"7LfM","chocolate":"d2691e","coral":"L7f50","cornflowerW":"6495ed","cornsilk":"Lf8dc","crimson":"dc143c","cyan":"MLL","NW":"MM8b","Ncyan":"M8b8b","Ngoldenrod":"b886b","Ngray":"a9a9a9","NV":"M64M","Nkhaki":"bdb76b","Nmagenta":"8bM8b","NoliveV":"556b2f","NY":"L8cM","Norchid":"9932cc","NS":"8bMM","Nsalmon":"e9967a","NseaV":"8fbc8f","NslateW":"483d8b","Nslategray":"2f4f4f","Nturquoise":"Mced1","Nviolet":"94Md3","OZ":"L1493","OskyW":"MbLf","dimgray":"696969","dodgerW":"1e90L","feldspar":"d19275","firebrick":"b22222","floralR":"Lfaf0","forestV":"228b22","fuchsia":"LML","gainsboro":"dcdcdc","ghostR":"f8f8L","gold":"Ld7M","goldenrod":"daa520","gray":"KKK","V":"M8M0","VX":"adL2f","honeydew":"f0Lf0","hotZ":"L69b4","indianS":"cd5c5c","indigo":"4bM82","ivory":"LLf0","khaki":"f0e68c","lavender":"e6e6fa","lavenderblush":"Lf0f5","lawnV":"7cfcM","lemonchiLon":"Lfacd","PW":"add8e6","Pcoral":"f0KK","Pcyan":"e0LL","PgoldenrodX":"fafad2","Pgrey":"d3d3d3","PV":"90K90","PZ":"Lb6c1","Psalmon":"La07a","PseaV":"20b2aa","PskyW":"87cefa","PslateW":"8470L","Pslategray":"778899","PstKlW":"b0c4de","PX":"LLe0","lime":"MLM","limeV":"32cd32","linen":"faf0e6","magenta":"LML","maroon":"8MM0","QUmarine":"66cdaa","QW":"MMcd","Qorchid":"ba55d3","Qpurple":"9370d8","QseaV":"3cb371","QslateW":"7b68K","QspringV":"Mfa9a","Qturquoise":"48d1cc","QvioletS":"c71585","midnightW":"191970","mintcream":"f5Lfa","mistyrose":"Le4e1","moccasin":"Le4b5","navajoR":"Ldead","navy":"MMK","oldlace":"fdf5e6","olive":"K8M0","olivedrab":"6b8e23","Y":"La5M","YS":"L45M","orchid":"da70d6","palegoldenrod":"Ke8aa","paleV":"98fb98","paleturquoise":"afKK","palevioletS":"d87093","papayawhip":"Lefd5","peachpuL":"Ldab9","peru":"cd853f","Z":"Lc0cb","plum":"dda0dd","powderW":"b0e0e6","purple":"8M0K","S":"LMM","rosybrown":"bc8f8f","royalW":"4169e1","saddlebrown":"8b4513","salmon":"faK72","sandybrown":"f4a460","seaV":"2e8b57","seashell":"Lf5K","sienna":"a0522d","silver":"c0c0c0","skyW":"87cKb","slateW":"6a5acd","slategray":"70K90","snow":"Lfafa","springV":"ML7f","stKlW":"4682b4","tan":"d2b48c","teal":"MKK","thistle":"d8bfd8","tomato":"L6347","turquoise":"40e0d0","violet":"K82K","violetS":"d02090","wheat":"f5deb3","R":"LLL","Rsmoke":"f5f5f5","X":"LLM","XV":"9acd32"}';
                        for (var i = 0; i < placeholder.length; i++) strEnColor = strEnColor.replace(new RegExp(placeholder[i], "g"), replace[i]);
                        return strEnColor;
                    },
                    __Hex2RGB = function (color) {
                        if (color.length === 4) {
                            var sColorNew = "#";
                            for (var i = 1; i < 4; i += 1) sColorNew += color.slice(i, i + 1).concat(color.slice(i, i + 1));
                            color = sColorNew;
                        }
                        var sColorChange = [];
                        for (var i = 1; i < 7; i += 2) sColorChange.push(parseInt("0x" + color.slice(i, i + 2)));
                        return "rgb(" + sColorChange[0] + "," + sColorChange[1] + "," + sColorChange[2] + ")";
                    };
                if (color.indexOf("rgba") >= 0) {
                    return color;
                } else if (color && /^#([0-9a-f]{3}|[0-9a-f]{6})$/.test(color)) {
                    return __Hex2RGB(color);
                } else if (color.indexOf("rgb") >= 0) {
                    if (color.indexOf("%") >= 0) {
                        var arrRgb = color.replace("rgb(", "").replace(")", "").replace(/%/g, "").split(",");
                        return "rgb(" + parseInt(arrRgb[0] / 100 * 255) + "," + parseInt(arrRgb[1] / 100 * 255) + "," + parseInt(arrRgb[2] / 100 * 255) + ")";
                    }
                    return color;
                } else {
                    var hexColor = JSON.parse(__EN_COLOR())[color];
                    if (!hexColor) {
                        console.error("抱歉！无该英文颜色");
                        return color;
                    }
                    return __Hex2RGB("#" + hexColor);
                }
            };
            var dimColor = function (color, opacity) {
                var opa = opacity || 0.5;
                if (color.indexOf("rgba") >= 0)
                    return color.substr(0, color.lastIndexOf(",")) + "," + parseFloat(color.substr(color.lastIndexOf(",") + 1, color.lastIndexOf(")") - 1)) * opa + ")";
                var arrRgb = toRGB(color).replace("rgb(", "").replace(")", "").split(",");
                return "rgba(" + arrRgb[0] + "," + arrRgb[1] + "," + arrRgb[2] + "," + opa + ")";
            };
            var isIE = function () {
                var userAgent = navigator.userAgent;
                var isIE = userAgent.indexOf("compatible") > -1 && userAgent.indexOf("MSIE") > -1; //判断是否IE<11浏览器
                var isEdge = userAgent.indexOf("Edge") > -1 && !isIE; //判断是否IE的Edge浏览器
                var isIE11 = userAgent.indexOf('Trident') > -1 && userAgent.indexOf("rv:11.0") > -1;
                return isIE || isIE11 || isEdge || false;
            };
            var isFire = function () {
                return navigator.userAgent.indexOf("Firefox") > -1;
            };
            var checkObjByType = function () {
                for (var i = 1; i < arguments.length; i++) if (typeof arguments[i] === arguments[0]) return arguments[i];
                return null;
            };
            var lzySelect = function (selector, copyStyle, callback) {
                this._randId = parseInt(Math.random() * 10000);
                this._hasNiceScroll = true;
                this._sOption = isFire() ? "option[selected=selected]" : "option:selected";
                this._size = 0;
                this._index = 0;
                this.distance = 2;
                this._callback = null;
                //-------判断参数，保证用户输入任何顺序都能允许-------------------------------------
                var realsSelector = checkObjByType("string", selector, copyStyle, callback);
                var isCopyStyle = checkObjByType("boolean", selector, copyStyle, callback);
                var realCallback = checkObjByType("function", selector, copyStyle, callback);
                //----初始化对象，-----------------------------------
                this._select = $(realsSelector);
                this._parent = $("<div class='lzy_select_parent'><div class='lzy_s_p_content'></div><div class='lzy_s_p_arraw'><svg><line x1='2.5' y1='12.5' x2='10' y2='7.5' style='stroke: rgb(0, 0, 0);'></line><line x1='10' y1='7.5' x2='17.5' y2='12.5' style='stroke: rgb(0, 0, 0);'></line></svg></div></div>")
                    .insertBefore(this._select.hide());
                var pStyles = ['position', 'top', 'right', 'left', 'bottom', 'float', 'margin-top', 'margin-left', 'margin-right', 'margin-bottom'];//位置方位的样式必须复制
                for (var i = 0; i < pStyles.length; i++) this._parent.css(pStyles[i], this._select.css(pStyles[i]));
                this._list = $("<div class='lzy_select_list'></div>").hide().appendTo("body");
                this._parent.find(".lzy_s_p_arraw").css({
                    "margin-top": this._parent.height() / 2 - this._parent.find(".lzy_s_p_arraw").height() / 2
                });
                var that = this;
                var themeColor = isCopyStyle ? that._select.css("border-left-color") : that._parent.css("border-left-color");

                if (isCopyStyle) {
                    //节点本身的样式根据是否添加来确定复制
                    var copyPStyle = ['width', 'height', 'color', 'backgroundColor', 'border-radius', 'border-left-color', 'border-right-color', 'border-top-color', 'border-bottom-color', 'border-style', 'border-width', 'box-shadow'];
                    for (i = 0; i < copyPStyle.length; i++) that._parent.css(copyPStyle[i], that._select.css(copyPStyle[i]));

                    that._parent.css({
                        'width': that._select.outerWidth() - that._select.css("border-width").replace("px", "") * 2,
                        'height': that._select.outerHeight() - that._select.css("border-width").replace("px", "") * 2
                    }).append("<style type=\"text/css\">.lzy_list_selected_" + that._randId + " > .selected {background-color: " + dimColor(themeColor, .12) + ";}" +
                        ".lzy_list_hover_" + that._randId + " > div:hover {background-color: " + dimColor(themeColor, .12) + ";}</style>")
                        .find("line").css({"stroke": themeColor});

                    var copyCStyle = ['paddingTop', 'paddingLeft', 'paddingRight', 'paddingBottom'];
                    for (i = 0; i < copyCStyle.length; i++) that._parent.find(".lzy_s_p_content").css(copyCStyle[i], that._select.css(copyCStyle[i]));
                    this._parent.find(".lzy_s_p_content").css({ //宽高要单独赋值是因为 css("width") 使用的是外宽，而实际要内宽
                        'width': that._select.width() - 23,
                        'height': that._select.height(),
                        'lineHeight': (isIE() ? that._select.height() + 4 : that._select.height()) + "px",
                        'color': that._select.css("color")
                    });

                    for (i = 3; i < copyPStyle.length; i++) that._list.css(copyPStyle[i], that._select.css(copyPStyle[i]));
                    that._list.addClass("lzy_list_hover_" + that._randId + " lzy_list_selected_" + that._randId).css('width', that._select.outerWidth() - that._select.css("border-width").replace("px", "") * 2);
                }

                var listSonWidth;
                try {
                    that._list.niceScroll({cursorcolor: themeColor, autohidemode: "hidden"});
                    listSonWidth = that._select.width() - 5 + "px";//5是 niceScroll 中滚动条的宽
                } catch (e) {
                    console.error("还没有导入 jquery-nicescroll.js！\n" + e);
                    that._hasNiceScroll = false;
                    that._parent.append("<style type=\"text/css\">.lzy_list_scroll_" + that._randId + "::-webkit-scrollbar-thumb{background-color:" + themeColor + ";}</style>");
                    that._list.addClass("lzy_list_scroll_" + that._randId);
                    listSonWidth = that._select.width() - (3 + 2) + "px";//3是 滚动条的宽，2是 边框的和
                    if (isIE() || isFire()) {
                        //火狐,ie 的滚动条的宽为15px
                        listSonWidth = that._select.width() - (16 + 2) + "px";//16是 滚动条的宽，2是 边框的和
                    }
                }

                var setIndex = function (index) {
                    that._list.find("div").eq(index).attr("selected", true).addClass("selected")
                        .siblings().removeClass("selected").removeAttr("selected");
                    var text = that._list.find("div").eq(index).text();
                    that._parent.attr("title", text).find(".lzy_s_p_content").html(text);
                };

                var setList = function () {
                    that._list.empty();
                    var options = that._select.find("option");
                    for (var i = 0; i < options.length; i++) {
                        var attrs = options[i].attributes;
                        var $son = $("<div></div>").attr("title", options[i].text).append(options[i].text);
                        if (isCopyStyle) {
                            $son.css({
                                'color': that._select.css("color"),
                                'padding-left': that._select.css("padding-left"),
                                'padding-right': that._select.css("padding-right"),
                                'padding-top': that._select.css("padding-top"),
                                'padding-bottom': that._select.css("padding-bottom"),
                                'line-height': that._select.height() + "px",
                                'height': that._select.height() + "px",
                                'width': listSonWidth
                            });
                        }
                        for (var j = 0; j < attrs.length; j++) {
                            if (attrs[j].name !== "class") $son.attr(attrs[j].name, attrs[j].value);
                        }
                        that._list.append($son);
                        if (i === 0) that._parent.attr("title", $son.attr("selected", true).text()).find(".lzy_s_p_content").html($son.attr("selected", true).text());
                    }
                    that._size = options.length;
                    setIndex(that._index);
                };
                setList();
                var setAnimate = function (type) {
                    if (type === "close") {
                        that._list.hide().removeClass("lzy_select_list_show");
                        if (that._hasNiceScroll) that._list.getNiceScroll().hide();
                    }
                    that._parent.find(".lzy_s_p_arraw").removeClass("lzy_s_p_arraw_" + (type === "close" ? "open" : "close")).addClass("lzy_s_p_arraw_" + type);
                };
                //----设置事件，-----------------------------------
                that._parent.on("click", function (event) {
                    (event || window.event).stopPropagation();
                    if (that._list.css("display") === "none") {
                        setAnimate("open");
                        var rect = that._parent.get(0).getBoundingClientRect();
                        var top = rect.bottom + that.distance;
                        if (that._list.height() > document.documentElement.clientHeight - rect.bottom) {//底部无足够空间
                            top = rect.top - that._list.outerHeight() - that.distance;
                        }
                        that._list.css({
                            "top": top,
                            "left": rect.left
                        }).show().addClass("lzy_select_list_show");
                        setTimeout(function () {
                            if (that._hasNiceScroll) that._list.getNiceScroll().show();
                        }, 200);
                    } else {
                        setAnimate("close");
                    }
                });

                var wheelTime;
                $(document).click(function () {
                    setAnimate("close");
                }).on("mousewheel DOMMouseScroll", function () {
                    //监听屏幕滚动时发生下拉宽移动响应
                    clearTimeout(wheelTime);
                    wheelTime = setTimeout(function () {
                        var rect = that._parent.get(0).getBoundingClientRect();
                        var top = that._list.height() > document.documentElement.clientHeight - rect.bottom ?
                            rect.top - that._list.outerHeight() - that.distance : rect.bottom + that.distance;
                        that._list.animate({"top": top}, "fast", "swing");
                    }, 200);
                });

                that._list.on("click", "div", function (event) {
                    (event || window.event).stopPropagation();
                    setAnimate("close");
                    that._parent.attr("title", $(this).attr("title")).find(".lzy_s_p_content").html($(this).text());
                    $(this).attr("selected", true).siblings().removeAttr("selected");
                    that._select.find("option").eq($(this).index()).attr("selected", true).siblings().attr("selected", false);
                    that._select.trigger("change");//引发change事件，保证原来的节点事件正常
                    if (realCallback) realCallback.call(that);
                    if (that._callback) that._callback.call(that);
                });

                /**
                 * 定时器，用于一直监听原 select 下拉框和 美化下拉框的内容是否一致
                 */
                setInterval(function () {
                    var size = that._select.find("option").length;
                    if (that._size !== size) setList();
                    var index = that._select.find(that._sOption).index();
                    if (that._index !== index) {
                        setIndex(index);
                        that._index = index;
                    }
                }, 50);
            };
            lzySelect.prototype = {
                change: function (callback) {
                    if (callback) this._callback = callback;
                },
                addOption: function (option) {
                    this._select.append(option);
                },
                addRow: function (text, params) {
                    if (!text || !params) throw new Error("传入的参数有误！");
                    var keys = Object.keys(params);
                    var $option = $("<option>" + text + "</option>");
                    for (var i = 0; i < keys.length; i++) {
                        $option.attr(keys[i], params[keys[i]]);
                    }
                    this._select.append($option);
                },
                getSelected: function (attrName) {
                    return this._select.find(this._sOption).attr(attrName || "value");
                },
                getSelectedText: function () {
                    return this._select.find(this._sOption).text();
                },
                getSelectedParams: function () {
                    var params = {};
                    var $option = this._select.find(this._sOption);
                    var attrs = $option.get(0).attributes;
                    for (var i = 0; i < attrs.length; i++) params[attrs[i].name] = attrs[i].value;
                    params.text = $option.text();
                    params.html = $option.text();
                    return params
                },
                setSelected: function (attrName, attrVal) {
                    if (!attrName && !attrVal) return false;
                    if (!attrVal) {
                        attrVal = attrName;
                        attrName = "value";
                    }
                    this._select.find("option[" + attrName + "='" + attrVal + "']").attr("selected", true).siblings().attr("selected", false);
                }
            };
            return new lzySelect(selector, copyStyle, callback);
        }
    });
})(jQuery, window, document);
