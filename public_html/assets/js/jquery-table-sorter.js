;(function ($, undefined) {

    var defaults = {
        arrows		: { up : '&#x25B4;', down : '&#x25BE;' },
		columns		: 'th',
		icon_class	: 'ts-sort-indicator'
    };

    function Plugin( element, options ) {
        this.element = $( element );
        this.options = $.extend( {}, defaults, options );
		this.isSorting = false;
        this.init();
    }

    Plugin.prototype.init = function () {
		var _self = this;
		this.element.find( this.options.columns ).on( 'click', function () {
			_self.onColumnClick.call( _self, this );
		});
    };
	
	Plugin.prototype.onColumnClick = function ( element ) {
		var $this = $( element ), index = $this.index();
		
		this.removeSortIndicators();
		this.setSortingIcon( element );
		this.sortColumn( index );
	};
	
	Plugin.prototype.removeSortIndicators = function () {
		this.element.find( '.' + this.options.icon_class ).remove();
	};
	
	Plugin.prototype.setSortingIcon = function ( column ) {
		$('<span />')
			.addClass( this.options.icon_class )
			.html((this.isSorting != true) ? this.options.arrows.up : this.options.arrows.down )
			.appendTo( column );
		
		this.isSorting = !this.isSorting;
		
		return this;
	};
	
	Plugin.prototype.sortColumn = function ( index ) {
		var _self = this, $trs = this.element.find( 'tbody tr' );
		
		$trs.sort( function ( tr1, tr2 ) {
			var v1 = _self.getTDValueByIndex( tr1, index ),
			    v2 = _self.getTDValueByIndex( tr2, index );
							
			return _self.compareValues.call( _self, v1, v2 );
		});
		
		this.element.find( 'tbody' ).empty().append( $trs );
		return this;
	};
	
	Plugin.prototype.compareValues = function ( v1, v2 ) {		
		if ( isNaN( v1 ) ) {
			if ( this.isSorting ) return v1.localeCompare( v2 );
			return v2.localeCompare( v1 );
		}
		else return this.sortNumbers( v1, v2 );
	};
	
	Plugin.prototype.sortNumbers = function ( v1, v2 ) {			
		var order  = ( v1 < v2 ) ? -1 : 0;
		
		if ( v1 > v2 ) order = 1;
		if ( this.isSorting && order !== 0 ) order = order * -1; 
		
		return order;
	};
	
	Plugin.prototype.getTDValueByIndex = function ( tr, index ) {
		var val = $.trim( $( tr ).find( 'td:eq( ' + index + ' )').text());
		return ( ! isNaN( val ) ) ? parseInt( val ) : val;
	};

    $.fn.TableSorter = function ( options ) {
        return this.each(function () {
            if ( ! $.data( this, "plugin" ) ) {
                $.data( this, "plugin", new Plugin( this, options ) );
            }
        });
    };
})( jQuery );