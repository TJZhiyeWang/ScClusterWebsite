function Bmove(oBtn,oImg){
		this.btn=$(oBtn).find("li");
		this.img=$(oImg).find("a");
		this.show();
	}
Bmove.prototype={
	show:function(){
			var thb=this;
			this.btn.each(function(i){			   
				$(this).mouseenter(function(){
					thb.btn.find("span").stop().animate({"width":"4px"});
					$(this).find("span").stop().animate({"width":"100%"});
					thb.img.eq(i).stop().animate({"opacity":"1","z-index":2});
					thb.img.eq(i).siblings("a").stop().animate({"opacity":"0","z-index":1});
					return false;
				})
			})	
		}
	}